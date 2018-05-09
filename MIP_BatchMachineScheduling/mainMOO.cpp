#include "General.h"
#include "Solver.h"

using namespace std;

int main ( int argc, char* argv[] ){

	//initialize problem instance from input file (demo if no input file available)
	char* inputFile;
	int eRates = 1;				// 1: Winter rates, 2: Summer rates
	int timeLimit = 3600;

	Problem problem;
	if (argc > 1) {
		problem.initializeFromFile(argv[1]);
		inputFile = argv[1];
		if (argc > 2) {
			eRates = atoi(argv[2]);
			if (argc > 3) {
				timeLimit = atoi(argv[3]);
			}
		}
	}
	else {
		problem.initializeFromFile("EEBS_demoProblem.dat");
		inputFile = "EESM_demoProblem.dat";
	}

	
	


	//EXPERIMENTS
	conductExperiments_WSC2018_specialCase(inputFile, problem, eRates, timeLimit);
	//conductExperiments_WSC2018_generalCase(inputFile, problem, eRates, timeLimit);


	
	//DEBUG
	//Solution test1 = solveTWCT_boundEPC(problem, 1, 258, 3600);
	//Solution test2 = solveEPC_boundTWCT(problem, 1, 555, 3600);	
	//Solution test3 = solveSimpleTWCT_boundEPC(problem, 1, 24000, 600);
	//Solution test4 = solveSimpleEPC_boundTWCT(problem, 1, 6666666, 600);
	
	return 0;
}




void generateInstances_WSC2018_specialCase() {
	int seed = 2052018;
	vector<int> numJobsPerF = {90, 150, 210}; 
	vector<int> batchSizes = {4, 8};
	vector<int> numFamilies = {3};
	vector<int> numMachines = {2, 3};
	vector<float> alphas = {0.0};					// 0 for r = 0
	vector<float> betas = { 0.5 };
	int replications = 5;

	Problem::generateInstances(seed, numJobsPerF, numFamilies, numMachines, batchSizes, alphas, betas, replications);

}
void generateInstances_WSC2018_specialCase_n50() {
	int seed = 2052018;
	vector<int> numJobsPerF = { 50 };
	vector<int> batchSizes = { 4, 8 };
	vector<int> numFamilies = { 3 };
	vector<int> numMachines = { 2, 3 };
	vector<float> alphas = { 0.0 };					// 0 for r = 0
	vector<float> betas = { 0.5 };
	int replications = 5;

	Problem::generateInstances(seed, numJobsPerF, numFamilies, numMachines, batchSizes, alphas, betas, replications);

}
void generateInstances_WSC2018_specialCase_n30() {
	int seed = 2052018;
	vector<int> numJobsPerF = { 30 };
	vector<int> batchSizes = { 4, 8 };
	vector<int> numFamilies = { 5 };
	vector<int> numMachines = { 2, 3 };
	vector<float> alphas = { 0.0 };					// 0 for r = 0
	vector<float> betas = { 0.5 };
	int replications = 5;

	Problem::generateInstances(seed, numJobsPerF, numFamilies, numMachines, batchSizes, alphas, betas, replications);

}
void generateInstances_WSC2018_specialCase_n35() {
	int seed = 2052018;
	vector<int> numJobsPerF = { 35 };
	vector<int> batchSizes = { 4, 8 };
	vector<int> numFamilies = { 3 };
	vector<int> numMachines = { 2, 3 };
	vector<float> alphas = { 0.0 };					// 0 for r = 0
	vector<float> betas = { 0.5 };
	int replications = 5;

	Problem::generateInstances(seed, numJobsPerF, numFamilies, numMachines, batchSizes, alphas, betas, replications);

}
void generateInstances_WSC2018_specialCase_n21() {
	int seed = 2052018;
	vector<int> numJobsPerF = { 21 };
	vector<int> batchSizes = { 4, 8 };
	vector<int> numFamilies = { 5 };
	vector<int> numMachines = { 2, 3 };
	vector<float> alphas = { 0.0 };					// 0 for r = 0
	vector<float> betas = { 0.5 };
	int replications = 5;

	Problem::generateInstances(seed, numJobsPerF, numFamilies, numMachines, batchSizes, alphas, betas, replications);

}
void generateInstances_WSC2018_specialCase_n30_F3() {
	int seed = 2052018;
	vector<int> numJobsPerF = { 30 };
	vector<int> batchSizes = { 4, 8 };
	vector<int> numFamilies = { 3 };
	vector<int> numMachines = { 2, 3 };
	vector<float> alphas = { 0.0 };					// 0 for r = 0
	vector<float> betas = { 0.5 };
	int replications = 5;

	Problem::generateInstances(seed, numJobsPerF, numFamilies, numMachines, batchSizes, alphas, betas, replications);

}
void generateInstances_WSC2018_specialCase_n18_F5() {
	int seed = 2052018;
	vector<int> numJobsPerF = { 18 };
	vector<int> batchSizes = { 4, 8 };
	vector<int> numFamilies = { 5 };
	vector<int> numMachines = { 2, 3 };
	vector<float> alphas = { 0.0 };					// 0 for r = 0
	vector<float> betas = { 0.5 };
	int replications = 5;

	Problem::generateInstances(seed, numJobsPerF, numFamilies, numMachines, batchSizes, alphas, betas, replications);

}
void generateInstances_WSC2018_specialCase_n15_F5() {
	int seed = 2052018;
	vector<int> numJobsPerF = { 15 };
	vector<int> batchSizes = { 4, 8 };
	vector<int> numFamilies = { 5 };
	vector<int> numMachines = { 2, 3 };
	vector<float> alphas = { 0.0 };					// 0 for r = 0
	vector<float> betas = { 0.5 };
	int replications = 5;

	Problem::generateInstances(seed, numJobsPerF, numFamilies, numMachines, batchSizes, alphas, betas, replications);
}
void generateInstances_WSC2018_specialCase_n25_F3() {
	int seed = 2052018;
	vector<int> numJobsPerF = { 25 };
	vector<int> batchSizes = { 4, 8 };
	vector<int> numFamilies = { 3 };
	vector<int> numMachines = { 2, 3 };
	vector<float> alphas = { 0.0 };					// 0 for r = 0
	vector<float> betas = { 0.5 };
	int replications = 5;

	Problem::generateInstances(seed, numJobsPerF, numFamilies, numMachines, batchSizes, alphas, betas, replications);
}
void generateInstances_WSC2018_generalCase() {
	int seed = 2052018;
	vector<int> numJobsPerF = { 5, 7 };
	vector<int> batchSizes = { 3, 4 };
	vector<int> numFamilies = { 2 };
	vector<int> numMachines = { 2 };
	vector<float> alphas = { 0.25, 0.5 };	
	vector<float> betas = { 0.5 };
	int replications = 5;

	Problem::generateInstances(seed, numJobsPerF, numFamilies, numMachines, batchSizes, alphas, betas, replications);

}
void generateInstances_WSC2018_generalCase_n5_F2() {
	int seed = 2052018;
	vector<int> numJobsPerF = { 5 };
	vector<int> batchSizes = { 3, 4 };
	vector<int> numFamilies = { 2 };
	vector<int> numMachines = { 2 };
	vector<float> alphas = { 0.25, 0.5 };
	vector<float> betas = { 0.5 };
	int replications = 5;

	Problem::generateInstances(seed, numJobsPerF, numFamilies, numMachines, batchSizes, alphas, betas, replications);

}
void conductExperiments_WSC2018_generalCase(char* filename, Problem problem, int eRates, int timeLimit) {
	SolutionSet solutions;
	solutions = eps_TWCT_EPC(problem, eRates, timeLimit);
	solutions.save(filename);
}
void conductExperiments_WSC2018_specialCase(char* filename, Problem problem, int eRates, int timeLimit) {
	SolutionSet solutions;
	solutions = eps_Simple_TWCT_EPC(problem, eRates, timeLimit);
	solutions.save(filename);
}