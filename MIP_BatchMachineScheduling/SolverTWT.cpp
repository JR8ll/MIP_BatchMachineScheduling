#include "Solver.h"


ILOSTLBEGIN

using namespace std;

typedef IloArray<IloBoolVarArray> IloBoolVarArray2;
typedef IloArray<IloBoolVarArray2> IloBoolVarArray3;

int solve_TWT(Problem problem, int ub_epc, int timeLimit) {

	//epsilon constraint parameters
	int N;					// N number of jobs
	int M;					// M number of machines
	int T;					// T number of time buckets
	int cap;				// batch capacity 
	int Q;					// penalty term (large integer)

	//set solving environment
	IloEnv env;

	//define parameters
	IloIntArray p(env);		// job processing times
	IloIntArray d(env);		// job due dates
	IloIntArray r(env);		// job release dates
	IloIntArray s(env);		// job sizes
	IloIntArray w(env);		// job weights
	IloIntArray f(env);		// job family

	IloIntArray pt(env);	// batch processing times
	IloIntArray st(env);	// batch earliest starting times
	IloIntArray fm(env);	// batch family

	IloIntArray e(env);		// electricity prices

	//set parameters
	//load job and batch parameters from problem instance
	M = problem.m;
	N = problem.p.size();

	cap = problem.k;
	Q = problem.q;
	for(int i = 0; i < N; i++) {
		p.add(problem.p[i]);
		d.add(problem.d[i]);
		r.add(problem.r[i]);
		s.add(problem.s[i]);
		w.add(problem.w[i]);
		f.add(problem.f[i]);
		pt.add(problem.p[i]);
		st.add(problem.r[i]);
		fm.add(problem.f[i]);
	}

	T = problem.e1.size();
	for(int i = 0; i < T; i++) {
		e.add(problem.e1[i]);
	}
	

	//set solving environment
	IloModel mod(env);

	// define + initialize decision variables
	IloBoolVarArray2 x(env);					// x[j][k]		job2batch assignment
	for(int i = 0; i < N; i++) {
		x.add(IloBoolVarArray(env, N));			// j = 1,..., N ; k = 1,..., F (=N)
	}

	IloArray <IloArray <IloArray<IloBoolVar> > > y (env, N);
	for(int k = 0; k < N; k++) {
		y[k] = IloArray<IloArray<IloBoolVar> > (env, M);
		for(int i = 0; i < M; i++) {
			y[k][i] = IloArray<IloBoolVar>  (env, T);
			for(int z = 0; z < T; z++) {
				y[k][i][z] = IloBoolVar(env);
			}
		}
	}

	IloNumVarArray c(env, N);					// c[j], j = 1,..., N			job completion times
	for(int i = 0; i < N; i++) {
		c[i] = IloNumVar(env);
	}
	IloNumVarArray ct(env, N);					// ct[k], k = 1,..., F (=N)		batch completion times
	for(int i = 0; i < N; i++) {
		ct[i] = IloNumVar(env);
	}
	

	//objective function
	 
	IloExpr twtObjective(env);
	IloExpr tardiness(env);
	for(int i = 0; i < N; i++) {
		tardiness = (IloMax(0, c[i] - d[i]));
		twtObjective += (tardiness * w[i]);
	}
	mod.add(IloMinimize(env, twtObjective));
	twtObjective.end();
	
	//alternative objective function: TWCT
	/*
	IloExpr twctObjective(env);
	for (int i = 0; i < N; i++) {
		twctObjective += c[i] * w[i];
	}
	mod.add(IloMinimize(env, twctObjective));
	twctObjective.end();
	*/
	//constraint sets: (2) EPC constraint
	IloExpr epc(env);
	for(int k = 0; k < N; k++){			//sum(k in Batches)
		for(int i = 0; i < M; i++){		//sum(i in Machines)
			for(int z = 0; z < T; z++){	//sum(t in Time)
				IloExpr ecost(env);
				int t_end = std::min(T, (int)z + (int)pt[k]);
				for(int h = z; h < t_end; h++) {
					ecost += e[h];	
				}
				epc += y[k][i][z] * ecost;
			}
		}
	}
	mod.add(epc <= ub_epc);
	epc.end();

	//constraint sets: (3) every job is assigned to one batch
	for(int j = 0; j < N; j++) {			//forall(j in Jobs)
		IloExpr con3(env);
		for(int k = 0; k < N; k++) {		//sum(k in Batches)
			con3 += x[j][k];
		}
		mod.add(con3 == 1);
		con3.end();
	}

	//constraint sets: (4) batches´ capacity constraint
	for(int i = 0; i < M; i++) {			//forall(i in Machines)
		for(int k = 0; k < N; k++) {		//forall(k in Batches)
			IloExpr con4(env);
			for(int j = 0; j < N; j++) {	//sum(j in Jobs)
				con4 += s[j] * x[j][k];
			}
			mod.add(con4 <= cap);
			con4.end();
		}
	}
	

	//constraint sets: (5) a job can only be assigned to a batch of the same family
	for(int j = 0; j < N; j++) {			//forall(j in Jobs)
		for(int k = 0; k < N; k++) {		//forall(k in Batches)
			IloExpr con5(env);
			con5 = x[j][k] * (fm[k] - f[j]);	
			mod.add(con5 == 0);
			con5.end();
		}
	}

	//constraint sets: (6) every batch is started no more than once
	for(int k = 0; k < N; k++) {					//forall(k in Batches)
		IloExpr con6(env);
		for(int i = 0; i < M; i++) {			//sum(i in Machines)
			for(int z = 0; z < T-1; z++) {	//sum(t in Time_1)
				con6 += y[k][i][z];
			}
		}
		mod.add(con6 <= 1);
		con6.end();
	}

	//constraint sets: (7) every batch is started at least once if it´s not empty
	for(int k = 0; k < N; k++) {					//forall(k in Batches)
		for(int j = 0; j < N; j++) {				//forall(j in Jobs)
			IloExpr con7(env);
			IloExpr sum7(env);
			for(int i = 0; i < M; i++) {			//sum(i in Machines)
				for(int z = 0; z < T-1; z++) {		//sum(t in Time_1)
					sum7 += y[k][i][z];
				}
			}
			con7 = sum7 - x[j][k];	
			mod.add(con7 >= 0);
			con7.end();
			sum7.end();
		}
	}

	//constraint sets: (8) every batch that is started is completed at or before the end of the shift
	for(int k = 0; k < N; k++) {					//forall(k in Batches)
		for(int j = 0; j < N; j++) {				//forall(j in Jobs)
			IloExpr con8(env);
			IloExpr sum8(env);
			for(int i = 0; i < M; i++) {			//sum(i in Machines)
				for(int z = 0; z < T-1; z++) {	//sum(t in Time_1)
					sum8 += (y[k][i][z] * ((IloNum)z + 1));			//ACHTUNG
				}
			}
			con8 = (-Q * (1 - x[j][k])) + sum8 + pt[k];
			mod.add(con8 <= T);
			con8.end();
			sum8.end();
		}
	}

	//constraint sets: (9) batches´ processing does not overlap on a machines
	for(int i = 0; i < M; i++) {					//forall(i in Machines)
		for(int z = 0; z < T; z++) {				//forall(t in Time)
			IloExpr con9(env);
			for(int k = 0; k < N; k++) {			//sum(k in Batches)
				int t_start = (std::max)(0, (z - problem.pt[k]+1));		//ACHTUNG
				for(int h = t_start; h <= z; h++) {	//sum(h in max(t-pt[k]+1, 1)..t)
					con9 += y[k][i][h];
				}
			}
			mod.add(con9 <= 1);
			con9.end();
		}
	}

	//constraint sets: (10) batch completion time
	for(int k = 0; k < N; k++) {				//forall(k in Batches)
		IloExpr con10(env);
		IloExpr sum10(env);
		for(int i = 0; i < M; i++) {			//sum(m in Machines)
			for(int z = 0; z < (T-1); z++) {	//sum(t in Time_1)
				sum10 += (y[k][i][z] * ((IloNum)z+1));	//equals the batches starting time	
			}
		}
		con10 = sum10 + pt[k] - ct[k];
		mod.add(con10 == 0);
		con10.end();
		sum10.end();
	}

	//constraint sets: (11) a batch may not be started before its earliest starting time
	for(int k = 0; k < N; k++) {							//forall(k in Batches)
		IloExpr con11(env);
		for(int i = 0; i < M; i++) {						//sum(i in Machines)
			for(int z = 0; z < (problem.st[k]-1); z++) {	//sum(t in 1..st[k]-1)		//ACHTUNG
				con11 += y[k][i][z];
			}
		}
		mod.add(con11 <= 0);
		con11.end();
	}

	//constraint sets: (12) a job can only be assigned to a batch if it´s started not before the job´s release time
	for(int j = 0; j < N; j++) {			//forall(j in Jobs)
		for(int k = 0; k < N; k++) {		//forall(k in Batches)
			IloExpr con12(env);
			con12 = ct[k] - pt[k] - (r[j] * x[j][k]);
			mod.add(con12 >= 0);
			con12.end();
		}
	}

	//constraint sets: (13) + (14) a job´s completion time is equal to its batch´s completion time
	for(int j = 0; j < N; j++) {			//forall(j in Jobs)
		for(int k = 0; k < N; k++) {		//forall(k in Batches)
			IloExpr con13(env);
			IloExpr con14(env);
			con13 = c[j] - ct[k] + (IloNum)Q * (1 - x[j][k]);
			con14 = c[j] - ct[k] - (IloNum)Q * (1 - x[j][k]);
			mod.add(con13 >= 0);
			mod.add(con14 <= 0);
			con13.end();
			con14.end();
		}
	}
	
	//constraint sets: (15) a job´s completion time is not negative
	for(int j = 0; j < N; j++) {		//forall(j in Jobs)
		IloExpr con15(env);
		con15 = c[j];
		mod.add(con15 >= 0);
	}

	

	
	//int stop; // DEBUG
	// set timelimit and solve
	IloCplex cplex(mod);
	cplex.setParam(IloCplex::TiLim, timeLimit);
	cplex.setOut(env.getNullStream());		//no screen output
	//cplex.exportModel("mip.lp");			//export model to file 

	try {
		cplex.solve();
		//cout << "CPLEX solving..." << endl;				//DEBUG
		}
	catch (IloException& e) {
		cerr << "Exception thrown by cplex.solve():";
		e.print(cerr);
		env.end();
		//cin >> stop;
		return 9999999;
		}

	// return the objective value
	try {
		cout << endl << "TWT: " << cplex.getObjValue() << endl;		//DEBUG
		cout << "rel. MIPgap (TWT): " << cplex.getMIPRelativeGap() << endl; //DEBUG
		
		
		for(int job = 0; job < N; job++) {
			IloNum wert = cplex.getValue(c[job]);
			//cout << "c[" << job+1 << "]= " << wert << ", ";	//DEBUG
		}
		//cout << endl;
		for(int batch = 0; batch < N; batch++) {
			IloNum wert = cplex.getValue(ct[batch]);
			//cout << "ct[" << batch+1 << "]= " << wert << ", "; //DEBUG
		}
		//cout << endl;
		for(int job = 0; job < N; job++) {
			for(int batch = 0; batch < N; batch++){
				IloNum wert = cplex.getValue(x[job][batch]);
				//cout << "x[" << job +1 << "][" << batch+1 << "] = " << (int)wert << ", "; //DEBUG
			}
			//cout << endl;
		}
		//cout << endl;
		for(int batch = 0; batch < N; batch++) {
			for(int machine = 0; machine < M; machine++) {
				for(int time = 0; time < T; time++) {
					IloNum wert = cplex.getValue(y[batch][machine][time]);
					//cout << "y[" << batch+1 << "][" << machine+1 << "][" << time+1 << "] = " << (int)wert << endl; //DEBUG
				}
			}
			//cout << endl;
			
		}


		//cin >> stop;
		int objective = (int) cplex.getObjValue();
		env.end();;
		return objective;
	}
	catch (IloException& e) {
		cerr << "No exact solution available." << endl;
		
		e.print(cout);
		//cin >> stop;
		env.end();
		return 999999;
	}
	env.end();

}

void epsilon_TWT_EPC(Problem problem, int timeLimit, char* file) {
	int epc = 0;
	int twt;
	int lb_epc;
	vector<pair<int, int>> nds;
	ofstream solutions;
	solutions.open("ParetoFront.csv", ios::app);
	solutions << file << endl;
	solutions << "TWT\t" << "EPC" << endl;
	//compute upper bound on EPC
	for (int t = 0; t < problem.e1.size(); t++) {
		epc += problem.e1[t];
	}
	epc = epc * problem.m;

	//compute lower bound on EPC: 330
	lb_epc = solve_EPC(problem, 999999, timeLimit);	//best EPC possible: 163

	//find initial solution
	twt = solve_TWT(problem, epc, timeLimit);		//best TWT possible: 396
	epc = solve_EPC(problem, twt, timeLimit);		//worst EPC being part of Pareto-front: 288
	solutions << twt << "\t" << epc << endl;
	nds.push_back(pair<int, int>(twt, epc));
	while (twt != 999999 && epc != 999999 && epc != lb_epc) {
		twt = solve_TWT(problem, epc - 1, timeLimit);	// twt(287) = 401
		epc = solve_EPC(problem, twt, timeLimit);		// epc(401) = 284
		solutions << twt << "\t" << epc << endl;
		nds.push_back(pair<int, int>(twt, epc));
	}
	cout << "epsilon-constraint procedure finished." << endl;

}