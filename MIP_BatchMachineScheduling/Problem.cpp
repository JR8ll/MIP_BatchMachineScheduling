#include "Problem.h"

int Problem::count = 0;
int Problem::seed = time(NULL);

Problem::Problem() {}
Problem::Problem(int numJobsPerFamily, int batchSize, int numFamilies, int numMachines, float alpha, float beta) {
	++Problem::count;
	this->n = numJobsPerFamily * numFamilies;
	this->m = numMachines;
	this->l = numFamilies;
	this->k = batchSize;
	this->rAlpha = alpha;
	this->dBeta = beta;

	int pMax = 0;
	int pSum = 0;
	int rMax = 0;
	vector<int> pFamily;

	for (unsigned i = 0; i < numFamilies; i++) {
		//set processing time per family randomly
		int pRandom = 1 + (rand() % 10);
		int myP = 0;
		if (pRandom <= 2) {
			myP = 2;
		}
		else if (pRandom <= 4) {
			myP = 4;
		}
		else if (pRandom <= 7) {
			myP = 10;
		}
		else if (pRandom <= 9) {
			myP = 16;
		}
		else if (pRandom <= 10) {
			myP = 20;
		}
		pFamily.push_back(myP);
		pSum += myP * numJobsPerFamily;
		if (myP > pMax) {
			pMax = myP;
		}
	}

	//for each job...
	for (unsigned i = 0; i < this->n; i++) {
		
		//set release date - pSum must be computed at this point
		int myR = 0;
		if (alpha > 0.01) {
			int ub_myR = (int)floor(alpha / this->k * pSum);
			myR = rand() % (ub_myR + 1);
			this->r.push_back(myR);
			if (myR > rMax) {
				rMax = myR;
			}
		}
		else {
			this->r.push_back(myR);
		}

		//set due date
		int myD = 0;
		if (beta > 0.01) {
			int ub_myD = myR + (int)floor(beta / this->k * pSum);
			myD = rand() % (ub_myD + 1);
			this->d.push_back(myD);
		}
		else {
			this->d.push_back(myD);
		}

		//set size
		this->s.push_back(1);

		//set weight randomly
		int wRandom = 1 + (rand() % 5);		// weight between 1 and 5
		this->w.push_back(wRandom);

		//set family
		int myF = (i / numJobsPerFamily) + 1;
		this->f.push_back(myF);		
		
		//set processing time per job dependent on the family
		this->p.push_back(pFamily[myF-1]);
	}

	//compute time horizon
	float tNumerator = 0;
	float tDenominator = 0;
	for (unsigned i = 0; i < numFamilies; i++) {
		tNumerator += ceil((float)numJobsPerFamily / (float) this->k) * (float) pFamily[i];
	}
	tDenominator = (float) this->m * (float) pMax;
	this->T = ceil(1.1 * (int)ceil(tNumerator / tDenominator)) * pMax + rMax;		// expand minimum time by 20%
	this->T = (this->T % 6 == 0) ? this->T : this->T + 6 -(this->T % 6);			// ceil to a multiple of six															

	//set e values for stepwise price functions
	int eOnPeak = 10;
	int ePartialPeak = 9;
	int eOffPeak = 8;
	
	//compute t values for stepwise price functions
	int tWinterOffPeak = this->T / 2;
	int tSummerPartialPeakOne = this->T / 3;
	int tSummerOffPeak = this->T / 2;
	int tSummerPartialPeakTwo = 5 * this->T / 6;

	//set second electricity price function (winter + summer rates)
	for (unsigned t = 0; t < this->T; t++) {
		//set first electricity price function (winter rates)
		if (t < tWinterOffPeak) {
			e1.push_back(eOnPeak);
		}
		else {
			e1.push_back(eOffPeak);
		}

		//set second electricity price function (summer rates)
		if (t < tSummerPartialPeakOne) {
			e2.push_back(eOnPeak);
		}
		else if (t < tSummerOffPeak) {
			e2.push_back(ePartialPeak);
		}
		else if (t < tSummerPartialPeakTwo) {
			e2.push_back(eOffPeak);
		}
		else {
			e2.push_back(ePartialPeak);
		}
	}



	
}

bool Problem::initializeFromFile(char* file) {
	ifstream input(file);
	if(!input) {
		cerr << file << " not found." << endl;
		return 0;
	}
	else {
		string dummy;
		//set parameters
		input >> dummy;			// "seed="
		input >> this->seed;
		input >> dummy;			// ";"
		input >> dummy;			// "N="
		input >> this->n;
		input >> dummy;			// ";"
		input >> dummy;			// "M="
		input >> this->m;
		input >> dummy;			// ";"
		input >> dummy;			// "F="
		input >> this->l;
		input >> dummy;			// ";"
		input >> dummy;			// "K="
		input >> this->k;
		input >> dummy;			// ";"
		input >> dummy;			// "T="
		input >> this->T;
		input >> dummy;			// ";"

		getline(input, dummy);	// empty line

		// processing times
		input >> dummy;			// "p=["
		int p;
		int pSum = 0;
		for(int j = 0; j < this->n; j++) {
			input >> p;
			pSum += p;
			this->p.push_back(p);
			this->pt.push_back(p);
		}
		input >> dummy;			// "];"

		// ready times
		input >> dummy;			// "r=["
		int r;
		int rMax = 0;
		for(int j = 0; j < this->n; j++) {
			input >> r;
			this->r.push_back(r);
			this->st.push_back(r);
			if(r > rMax) rMax = r;
		}
		input >> dummy;			// "];"

		// due dates
		input >> dummy;			// "d=["
		int d;
		for(int j = 0; j < this->n; j++) {
			input >> d;
			this->d.push_back(d);
		}
		input >> dummy;			// "];"

		// sizes
		input >> dummy;			// "s=["
		int s;
		int sSum = 0;
		for(int j = 0; j < this->n; j++) {
			input >> s;
			sSum += s;
			this->s.push_back(s);
		}
		input >> dummy;			// "];"

		// weights
		input >> dummy;			// "w=["
		int w;
		for(int j = 0; j < this->n; j++) {
			input >> w;
			this->w.push_back(w);
		}
		input >> dummy;			// "];"

		// families
		input >> dummy;			// "f=["
		int f;
		for(int j = 0; j < this->n; j++) {
			input >> f;
			this->f.push_back(f);
		}
		input >> dummy;			// "];"

		getline(input, dummy);	// empty line

		// 1st electricity price
		input >> dummy;			// "e1=["
		
		int e;
		for(int t = 0; t < this->T; t++) {
			input >> e;
			this->e1.push_back(e);
		}
		input >> dummy;			// "];"

		// 2nd electricity price
		input >> dummy;			// "e2=["
		for (int t = 0; t < this->T; t++) {
			input >> e;
			this->e2.push_back(e);
		}
		input >> dummy;			// "];"

		this->q = pSum + rMax + this->T;	

		cout << "Problem instance from file " << file << " initialized." << endl << endl;
		return 1;
	};	
}
void Problem::saveToFile(int replication) {

	string path = "pInstances//";
	string file = "EEBS";
	string nString = to_string(this->n / this->l);
	string lString = to_string(this->l);
	string mString = to_string(this->m);
	string kString = to_string(this->k);
	int alph = (int)(this->rAlpha * 100.0);
	string aString = to_string(alph);
	int bet = (int)(this->dBeta * 100.0);
	string bString = to_string(bet);
	string repString = to_string(replication);
	string ul = to_string('_');
	string hy = to_string('-');

	ostringstream os;
	os << file << "_n" << nString << "_F" << lString << "_M" << mString << "_B" << kString << "_a0." << aString << "_b0." << bString << "T1.1" << "-" << repString << ".dat";
	string filename = os.str();
	ofstream out;
	out.open(path + filename.c_str());

	out << "seed= " << Problem::seed << " ;" << endl;
	out << "N= " << this->n << " ;" << endl;
	out << "M= " << this->m << " ;" << endl;
	out << "F= " << this->l << " ;" << endl;
	out << "K= " << this->k << " ;" << endl;
	out << "T= " << this->T << " ;" << endl;
	out << endl;
	out << "p=[";
	for (unsigned i = 0; i < this->n; i++) {
		out << " " << this->p[i];
	}
	out << " ];" << endl;
	out << "r=[";
	for (unsigned i = 0; i < this->n; i++) {
		out << " " << this->r[i];
	}
	out << " ];" << endl;
	out << "d=[";
	for (unsigned i = 0; i < this->n; i++) {
		out << " " << this->d[i];
	}
	out << " ];" << endl;
	out << "s=[";
	for (unsigned i = 0; i < this->n; i++) {
		out << " " << this->s[i];
	}
	out << " ];" << endl;
	out << "w=[";
	for (unsigned i = 0; i < this->n; i++) {
		out << " " << this->w[i];
	}
	out << " ];" << endl;
	out << "f=[";
	for (unsigned i = 0; i < this->n; i++) {
		out << " " << this->f[i];
	}
	out << " ];" << endl;
	out << endl;
	out << "e1=[";
	for (unsigned t = 0; t < this->T; t++) {
		out << " " << e1[t];
	}
	out << " ];" << endl;
	out << "e2=[";
	for (unsigned t = 0; t < this->T; t++) {
		out << " " << e2[t];
	}
	out << " ];";
	out.close();
}
void Problem::generateInstances(int in_seed, vector<int> nPerF, vector<int> numFamilies, vector<int> numMachines, vector<int> batchSizes, vector<float> alphas, vector<float> betas, int replications) {
	if (in_seed > 0) {
		Problem::seed = in_seed;
	}
	srand(Problem::seed);

	for (unsigned i = 0; i < nPerF.size(); i++) {
		for (unsigned j = 0; j < numFamilies.size(); j++) {
			for (unsigned k = 0; k < numMachines.size(); k++) {
				for (unsigned l = 0; l < batchSizes.size(); l++) {
					for (unsigned m = 0; m < alphas.size(); m++) {
						for (unsigned q = 0; q < betas.size(); q++) {
							for (unsigned n = 0; n < replications; n++) {
								Problem p = Problem(nPerF[i], batchSizes[l], numFamilies[j], numMachines[k], alphas[m], betas[q]);
								p.saveToFile(n + 1);
							}
						}
					}
				}
			}
		}
	}
}

bool Problem::hasUniformS() {
	if (this->s.empty()) {
		return false;	// return false if there are no jobs
	}
	else {
		int firstS = this->s[0];
		for (unsigned i = 1; i < this->s.size(); i++) {
			if (this->s[i] != firstS) {
				return false;  // return false if there are two jobs with different sizes
			}
		}
		return true;	// return true if there are jobs, and all sizes are equal
	}
}
bool Problem::hasUniformW() {
	if (this->w.empty()) {
		return false;	// return false if there are no jobs
	}
	else {
		int firstW = this->w[0];
		for (unsigned i = 1; i < this->w.size(); i++) {
			if (this->w[i] != firstW) {
				return false;  // return false if there are two jobs with different weights
			}
		}
		return true;	// return true if there are jobs, and all weights are equal
	}
}
bool Problem::hasUniformP() {
	if (this->p.empty()) {
		return false;	// return false if there are no jobs
	}
	else {
		int firstP = this->p[0];
		for (unsigned i = 1; i < this->p.size(); i++) {
			if (this->p[i] != firstP) {
				return false;  // return false if there are two jobs with different processing times
			}
		}
		return true;	// return true if there are jobs, and all processing times are equal
	}
}
bool Problem::hasSingleF() {
	if (this->f.empty()) {
		return false;	// return false if there are no jobs
	}
	else {
		int firstF = this->f[0];
		for (unsigned i = 1; i < this->f.size(); i++) {
			if (this->f[i] != firstF) {
				return false;  // return false if there are two jobs with a different family
			}
		}
		return true;	// return true if there are jobs, and all belong to the same family
	}
}
bool Problem::hasUniformR() {
	if (this->r.empty()) {
		return false;	// return false if there are no jobs
	}
	else {
		int firstR = this->r[0];
		for (unsigned i = 1; i < this->r.size(); i++) {
			if (this->r[i] != firstR) {
				return false;  // return false if there are two jobs with a different family
			}
		}
		return true;	// return true if there are jobs, and all belong to the same family
	}
}