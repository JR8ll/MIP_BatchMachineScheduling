#pragma once

#define NOMINMAX
#include <ilcplex/ilocplex.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <chrono>
#include <string>

//general case instances, Winter Simulation Conference 2018	
//void generateInstances_WSC2018_generalCase();		// declaration moved to solver.h		
//void conductExperiments_WSC2018_generalCase();

//special case instances (r=0), Winter Simulation Conference 2018
//void generateInstances_WSC2018_specialCase();		 // declaration moved to solver.h
//void conductExperiments_WSC2018_specialCase();

bool compDesc(int& a, int& b);