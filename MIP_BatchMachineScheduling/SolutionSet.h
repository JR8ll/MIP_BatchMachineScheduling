#pragma once
#include "General.h"

using namespace std;

struct Solution {
	string objective;
	bool solved;
	int value;
	int bound;
	double time;
	double relMIPgap;
	int eRates;
};

class SolutionSet {
public:
	vector< pair<Solution, Solution > > solutions;

	void save(char* filename);
};
