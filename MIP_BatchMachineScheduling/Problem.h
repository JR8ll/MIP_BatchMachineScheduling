#include "General.h"

using namespace std;

class Problem {

public: 
	static int count;	
	static int seed;

	int n;				// number of jobs
	int m;				// number of machines
	int l;				// number of families
	int k;				// batch capacity
	int T;				// time horizon
	int q;				// penalty term (large integer)

	float rAlpha;		// alpha factor for computation of r
	float dBeta;		// beta factor for computation of d
	
	vector<int> p;		// job processing times	
	vector<int> d;		// job due dates
	vector<int> r;		// job release dates
	vector<int> s;		// job sizes
	vector<int> w;		// job weights
	vector<int> f;		// job family

	vector<int> pt;		// batch processing times
	vector<int> st;		// batch earliest starting times

	vector<int> e1;		// electricity prices (winter tariff)
	vector<int> e2;		// electricity prices (summer tariff)
	
	static void generateInstances(int seed, vector<int> nPerF, vector<int> numFamilies, vector<int> numMachines, vector<int> batchSizes, vector<float> alphas, vector<float> betas, int replications);

	Problem();
	Problem(int numJobsPerFamily, int batchSize, int numFamilies, int numMachines, float alpha, float beta);

	bool initializeFromFile(char* file);
	void saveToFile(int replication);

	bool hasUniformS();
	bool hasUniformW();
	bool hasUniformP();
	bool hasSingleF();
	bool hasUniformR();


};