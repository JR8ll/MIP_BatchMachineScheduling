#include "Problem.h"
#include "SolutionSet.h"

int solve_TWT(Problem problem, int ub_epc, int timeLimit);						// SolverTWT.cpp 
int solve_EPC(Problem problem, int ub_twt, int timeLimit);						// SolverEPC.cpp
int solve_TwCT(Problem problem, int ub_epc, int eRates, int timeLimit);			// SolverTwCT.cpp TODO: verify/test
int solve_TCT(Problem problem, int ub_epc, int eRates, int timeLimit);			// SolverTwCT.cpp TODO: verify/test
void epsilon_TWT_EPC(Problem problem, int timeLimit, char* file);				// SolverTWT.cpp
void epsilon_TwCT_EPC(Problem problem, int timeLimit, char* file);				// SolverTwCT.cpp TODO: implement

Solution solveTWCT_boundEPC(Problem problem, int eRates, int ub_epc, int timeLimit);
Solution solveEPC_boundTWCT(Problem problem, int eRates, int ub_twct, int timeLimit);
SolutionSet eps_TWCT_EPC(Problem problem, int eRates, int timeLimit);			// SolverTwCT.cpp

Solution solveSimpleTWCT_boundEPC(Problem problem, int eRates, int ub_epc, int timeLimit);
Solution solveSimpleEPC_boundTWCT(Problem problem, int eRates, int ub_twct, int timeLimit);
SolutionSet eps_Simple_TWCT_EPC(Problem problem, int eRates, int timeLimit);

//general case instances, Winter Simulation Conference 2018	
void generateInstances_WSC2018_generalCase();
void generateInstances_WSC2018_generalCase_n5_F2();
void conductExperiments_WSC2018_generalCase(char* filename,  Problem problem, int eRates, int timeLimit);

//special case instances (r=0), Winter Simulation Conference 2018
void generateInstances_WSC2018_specialCase();
void generateInstances_WSC2018_specialCase_n50();
void generateInstances_WSC2018_specialCase_n30();
void generateInstances_WSC2018_specialCase_n35();
void generateInstances_WSC2018_specialCase_n21();
void generateInstances_WSC2018_specialCase_n18_F5();	// N=90
void generateInstances_WSC2018_specialCase_n30_F3();	// N=90
void generateInstances_WSC2018_specialCase_n15_F5();	// N=75
void generateInstances_WSC2018_specialCase_n25_F3();	// N=75
void conductExperiments_WSC2018_specialCase(char* filename, Problem problem, int eRates, int timeLimit);


// TODO constraint sets (11), ggf. weitere anpassen, EPC Versionen mit ub_tct, ub_twct implementieren, alternative eRates berücksichtigen, epsilons