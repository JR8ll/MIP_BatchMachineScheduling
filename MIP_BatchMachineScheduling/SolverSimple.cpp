#include "Solver.h"

ILOSTLBEGIN

using namespace std;

typedef IloArray<IloBoolVarArray> IloBoolVarArray2;
typedef IloArray<IloBoolVarArray2> IloBoolVarArray3;

Solution solveSimpleTWCT_boundEPC(Problem problem, int eRates, int ub_epc, int timeLimit) {
	//epsilon constraint parameters
	IloInt B;					// N number of batches
	IloInt M;					// M number of machines
	IloInt T;					// T number of time buckets

	//set solving environment
	IloEnv env;

	//set solving environment
	IloModel mod(env);

	//define parameters
	IloIntArray p(env);			// batch processing times
	IloIntArray w(env);			// batch weights (sum of assigned jobs´ weights)
	IloIntArray e1(env);		// electricity prices winter rates
	IloIntArray e2(env);		// electricity prices summer rates

	//form batches
	//TODO form batches

	vector< vector < int > > wBatches(problem.l);
	vector<int> pFam(problem.l);

	for (unsigned i = 0; i < problem.n; i++) {
		wBatches[problem.f[i] - 1].push_back(problem.w[i]);
		pFam[problem.f[i] - 1] = problem.p[i];
	}

	for (unsigned i = 0; i < problem.l; i++) {
		sort(wBatches[i].begin(), wBatches[i].end(), compDesc);
	}

	for (unsigned i = 0; i < problem.l; i++) {
		IloInt sumWeight = 0;
		for (unsigned j = 0; j < wBatches[i].size(); j++) {
			sumWeight += wBatches[i][j];
			if ((j + problem.k + 1) % problem.k == 0 || j == wBatches[i].size() - 1) {
				w.add(sumWeight);
				p.add(pFam[i]);
				sumWeight = 0;
			}
		}
	}


	//prepare solution container
	Solution sol;
	sol.objective = "TwCT";
	sol.bound = ub_epc;
	sol.eRates = eRates;

	//initialize parameters
	B = (IloInt) w.getSize();			
	M = (IloInt)problem.m;
	T = (IloInt)problem.T;

	for (unsigned t = 0; t < problem.T; t++) {
		e1.add((IloInt)problem.e1[t]);
		e2.add((IloInt)problem.e2[t]);
	}

	//DEBUG
	/*
	cout << "B (# of batches): " << B << endl;
	cout << "{";
	for (unsigned i = 0; i < B; i++) {
		cout << "(p" << p[i] << ", w" << w[i] << ")";
	}
	cout << "}" << endl;
	*/

	// define + initialize decision variables
	IloArray <IloArray <IloArray<IloBoolVar> > > x(env, B);
	for (int b = 0; b < B; b++) {
		x[b] = IloArray<IloArray<IloBoolVar> >(env, T);
		for (int t = 0; t < T; t++) {
			x[b][t] = IloArray<IloBoolVar>(env, M);
			for (int m = 0; m < M; m++) {
				x[b][t][m] = IloBoolVar(env);
			}
		}
	}

	//objective function: TwCT
	IloExpr twctObjective(env);
	for (int b = 0; b < B; b++) {
		for (int t = 0; t < T; t++) {
			for (int m = 0; m < M; m++) {
				twctObjective += ((IloInt)t + p[b]) * w[b] * x[b][t][m];
			}
		}
	}
	mod.add(IloMinimize(env, twctObjective));

	//constraint set (1): EPC constraint
	IloExpr epc1(env);
	IloExpr epc2(env);
	for (int b = 0; b < B; b++) {			// sum(b in Batches)
		for (int t = 0; t < T; t++) {		// sum(t in Time)
			for (int m = 0; m < M; m++) {	// sum(m in Machines)
				IloExpr ecost1(env);
				IloExpr ecost2(env);
				for (int h = t; h < min((int)T, t + (int)p[b]); h++) {
					ecost1 += e1[h];
					ecost2 += e2[h];
				}
				epc1 += x[b][t][m] * ecost1;
				epc2 += x[b][t][m] * ecost2;
			}
		}
	}
	if (eRates == 1) {
		mod.add(epc1 <= (IloInt)ub_epc);
	}
	else {
		mod.add(epc2 <= (IloInt)ub_epc);
	}
	epc1.end();
	epc2.end();
	
	//constraint set (2): every batch is processed once
	for (int b = 0; b < B; b++) {					// forall(b in Batches)
		IloExpr con2(env);	
		for (int t = 0; t < (T - p[b]); t++) {		// sum(t in 0..T-pt[b])
			for (int m = 0; m < M; m++) {			// sum(m in Machines)
				con2 += x[b][t][m];
			}
		}
		mod.add(con2 == 1);
	}

	//constraint set (3): non-overlapping batch processing
	for (int m = 0; m < M; m++) {										// forall(m in Machines)
		for (int t = 0; t < T; t++) {									// forall(t in Time)
			IloExpr con3(env);
			for (int b = 0; b < B; b++) {								// sum(b in Batches)
				for (int h = std::max(t -(int)p[b], 0); h < t; h++) {		// sum(h in t-p[b]..T)
					con3 += x[b][h][m];
				}
			}
			mod.add(con3 <= 1);
		}
	}

	//constraint set (4): batches are to be completed before T
	for (int b = 0; b < B; b++) {										// forall(b in Batches)
		IloExpr con4(env);
		for (int t = T - p[b]; t < T; t++) {							// sum(t in (T-p[b])..T)
			for (int m = 0; m < M; m++) {								// sum(m in Machines)
				con4 += x[b][t][m];
			}
		}
		mod.add(con4 <= 0);
	}

	//set timelimit and solve
	IloCplex cplex(mod);
	cplex.setParam(IloCplex::TiLim, timeLimit);
	cplex.setOut(env.getNullStream());

	try {
		IloNum start = cplex.getCplexTime();
		cplex.solve();
		sol.time = cplex.getCplexTime() - start;
	}
	catch (IloException& e) {
		cerr << "Exception thrown by cplex.solve():";
		e.print(cerr);
		env.end();
		sol.solved = false;
		return sol;
	}

	// return the objective value
	try {

		//DEBUGGING information for cout
		/*
		cout << endl << "TwCT: " << (int) round(cplex.getObjValue()) << endl;		//DEBUG
		cout << "rel. MIPgap (TwCT): " << cplex.getMIPRelativeGap() << endl; //DEBUG

		cout << endl;
		for (int batch = 0; batch < B; batch++) {
		for (int time = 0; time < T; time++) {
		for (int machine = 0; machine < M; machine++) {
		IloInt wert = cplex.getValue(x[batch][time][machine]);
		cout << "x[" << batch+1 << "][" << time << "][" << machine+1 << "] = " << (int)wert << endl; //DEBUG
		}
		}
		cout << endl;

		}
		*/

		sol.value = (int)round(cplex.getObjValue());
		sol.relMIPgap = (float)cplex.getMIPRelativeGap();
		sol.solved = true;
		env.end();
		return sol;
	}
	catch (IloException& e) {
		cerr << "No exact solution available." << endl;
		e.print(cout);
		env.end();
		sol.solved = false;
		return sol;
	}
	env.end();
}

Solution solveSimpleEPC_boundTWCT(Problem problem, int eRates, int ub_twct, int timeLimit) {
	
	//define parameters
	IloInt B;					// N number of batches
	IloInt M;					// M number of machines
	IloInt T;					// T number of time buckets
	
	//set solving environment
	IloEnv env;								
	
	//set solving environment
	IloModel mod(env);

	IloIntArray p(env);			// batch processing times
	IloIntArray w(env);			// batch weights (sum of assigned jobs´ weights)
	IloIntArray e1(env);		// electricity prices winter rates
	IloIntArray e2(env);		// electricity prices summer rates

	//form batches
	//TODO form batches
	
	vector< vector < int > > wBatches(problem.l);
	vector<int> pFam(problem.l);

	for (unsigned i = 0; i < problem.n; i++) {
		wBatches[problem.f[i] - 1].push_back(problem.w[i]);
		pFam[problem.f[i] - 1] = problem.p[i];
	}

	for (unsigned i = 0; i < problem.l; i++) {
		sort(wBatches[i].begin(), wBatches[i].end(), compDesc);
	}

	for (unsigned i = 0; i < problem.l; i++) {
		IloInt sumWeight = 0;
		for (unsigned j = 0; j < wBatches[i].size(); j++) {
			sumWeight += wBatches[i][j];
			if ((j + problem.k + 1) % problem.k == 0 || j == wBatches[i].size() - 1) {
				w.add(sumWeight);
				p.add(pFam[i]);
				sumWeight = 0;
			}
		}
	}

	//prepare solution container
	Solution sol;
	sol.objective = "EPC";
	sol.bound = ub_twct;
	sol.eRates = eRates;

	//initialize parameters
	B = (IloInt) w.getSize();			// TODO compute from batch formation
	M = (IloInt)problem.m;
	T = (IloInt)problem.T;

	for (unsigned t = 0; t < problem.T; t++) {
		e1.add((IloInt)problem.e1[t]);
		e2.add((IloInt)problem.e2[t]);
	}

	// define + initialize decision variables
	IloArray <IloArray <IloArray<IloBoolVar> > > x(env, B);
	for (int b = 0; b < B; b++) {
		x[b] = IloArray<IloArray<IloBoolVar> >(env, T);
		for (int t = 0; t < T; t++) {
			x[b][t] = IloArray<IloBoolVar>(env, M);
			for (int m = 0; m < M; m++) {
				x[b][t][m] = IloBoolVar(env);
			}
		}
	}

	//objective function: EPC objective
	IloExpr epc1(env);
	IloExpr epc2(env);
	for (int b = 0; b < B; b++) {			// sum(b in Batches)
		for (int t = 0; t < T; t++) {		// sum(t in Time)
			for (int m = 0; m < M; m++) {	// sum(m in Machines)
				IloExpr ecost1(env);
				IloExpr ecost2(env);
				for (int h = t; h < min((int)T, t + (int)p[b]); h++) {
					ecost1 += e1[h];
					ecost2 += e2[h];
				}
				epc1 += x[b][t][m] * ecost1;
				epc2 += x[b][t][m] * ecost2;
			}
		}
	}
	if (eRates == 1) {
		mod.add(IloMinimize(env, epc1));
	}
	else {
		mod.add(IloMinimize(env, epc2));
	}
	epc1.end();
	epc2.end();

	//constraint set (1): TwCT
	IloExpr twctConstraint(env);
	for (int b = 0; b < B; b++) {
		for (int t = 0; t < T; t++) {
			for (int m = 0; m < M; m++) {
				twctConstraint += ((IloInt)t + p[b]) * w[b] * x[b][t][m];
			}
		}
	}
	mod.add(twctConstraint <= ub_twct);

	//constraint set (2): every batch is processed once
	for (int b = 0; b < B; b++) {					// forall(b in Batches)
		IloExpr con2(env);
		for (int t = 0; t < (T - p[b]); t++) {		// sum(t in 0..T-pt[b])
			for (int m = 0; m < M; m++) {			// sum(m in Machines)
				con2 += x[b][t][m];
			}
		}
		mod.add(con2 == 1);
	}

	//constraint set (3): non-overlapping batch processing
	for (int m = 0; m < M; m++) {										// forall(m in Machines)
		for (int t = 0; t < T; t++) {									// forall(t in Time)
			IloExpr con3(env);
			for (int b = 0; b < B; b++) {								// sum(b in Batches)
				for (int h = std::max(t - (int)p[b], 0); h < t; h++) {		// sum(h in t-p[b]..T)
					con3 += x[b][h][m];
				}
			}
			mod.add(con3 <= 1);
		}
	}

	//constraint set (4): batches are to be completed before T
	for (int b = 0; b < B; b++) {										// forall(b in Batches)
		IloExpr con4(env);
		for (int t = T - p[b]; t < T; t++) {							// sum(t in (T-p[b])..T)
			for (int m = 0; m < M; m++) {								// sum(m in Machines)
				con4 += x[b][t][m];
			}
		}
		mod.add(con4 <= 0);
	}

	//set timelimit and solve
	IloCplex cplex(mod);
	cplex.setParam(IloCplex::TiLim, timeLimit);
	cplex.setOut(env.getNullStream());

	try {
		IloNum start = cplex.getCplexTime();
		cplex.solve();
		sol.time = cplex.getCplexTime() - start;
	}
	catch (IloException& e) {
		cerr << "Exception thrown by cplex.solve():";
		e.print(cerr);
		env.end();
		sol.solved = false;
		return sol;
	}

	// return the objective value
	try {

		//DEBUGGING information for cout
		/*
		cout << endl << "TwCT: " << (int) round(cplex.getObjValue()) << endl;		//DEBUG
		cout << "rel. MIPgap (TwCT): " << cplex.getMIPRelativeGap() << endl; //DEBUG

		cout << endl;
		for (int batch = 0; batch < B; batch++) {
		for (int time = 0; time < T; time++) {
		for (int machine = 0; machine < M; machine++) {
		IloInt wert = cplex.getValue(x[batch][time][machine]);
		cout << "x[" << batch+1 << "][" << time << "][" << machine+1 << "] = " << (int)wert << endl; //DEBUG
		}
		}
		cout << endl;

		}
		*/

		sol.value = (int)round(cplex.getObjValue());
		sol.relMIPgap = (float)cplex.getMIPRelativeGap();
		sol.solved = true;
		env.end();
		return sol;
	}
	catch (IloException& e) {
		cerr << "No exact solution available." << endl;
		e.print(cout);
		env.end();
		sol.solved = false;
		return sol;
	}
	env.end();
}

SolutionSet eps_Simple_TWCT_EPC(Problem problem, int eRates, int timeLimit) {

	if (!problem.hasUniformS() || !problem.hasUniformR()) {
		cout << "+++ WARNING: Problem is not applicable for simplified MIP +++" << endl;
	}
	SolutionSet sols;
	int ub_epc;
	int ub_epc1 = 0;
	int ub_epc2 = 0;
	int solCount = 0;
	for (unsigned t = 0; t < problem.T; t++) {
		ub_epc1 += problem.e1[t];
		ub_epc2 += problem.e2[t];
	}
	ub_epc1 = ub_epc1 * problem.m;
	ub_epc2 = ub_epc2 * problem.m;
	if (eRates == 1) {
		ub_epc = ub_epc1;
	}
	else {
		ub_epc = ub_epc2;
	}

	pair<Solution, Solution> point;

	do {
		Solution sol1 = solveSimpleTWCT_boundEPC(problem, eRates, ub_epc, timeLimit);
		Solution sol2 = solveSimpleEPC_boundTWCT(problem, eRates, sol1.value, timeLimit);
		if (sol1.solved == false || sol2.solved == false) {
			cout << "Solution set completed." << endl;
			break;
		}
		else {
			solCount++;
			point.first = sol1;
			point.second = sol2;
			sols.solutions.push_back(point);
			cout << "Solution " << solCount << " added (" << sol1.objective << "=" << sol1.value << "[b" << sol1.bound << "] | " << sol2.objective << "=" << sol2.value << "[b" << sol2.bound << "])." << endl;
			ub_epc = sol2.value;
			ub_epc--;
		}
	} while (true);		// Exit through the break statement

	return sols;
}


bool compDesc(int& a, int& b) {
	return a > b;
}