#include "SolutionSet.h"

void SolutionSet::save(char* filename) {
	
	// open file
	string file = "EEBS_solutions.csv";
	ofstream out;
	out.open(file.c_str(), ios::app);

	if (solutions.empty()) {
		out << "Problem\t" << filename << "\t" << "no solutions" << endl;
	}
	else {
		// sum up computation times
		double sumCompTime = 0.0;
		for (unsigned i = 0; i < solutions.size(); i++) {
			sumCompTime += solutions[i].first.time + solutions[i].second.time;
		}

		// write summary
		out << "Problem\t" << filename << "\t" << solutions.size() << "\t" << solutions[0].first.eRates << "\t" << sumCompTime << endl;

		// write headings
		out << solutions[0].first.objective << "\t" << solutions[0].second.objective << "\t" << "t_" << solutions[0].first.objective << "\t" << "t_" << solutions[0].second.objective << "\t" << "gap_" << solutions[0].first.objective << "\t" << "gap_" << solutions[0].second.objective << endl;

		// write solutions
		for (unsigned i = 0; i < solutions.size(); i++) {
			out << solutions[i].first.value << "\t" << solutions[i].second.value << "\t" << solutions[i].first.time << "\t" << solutions[i].second.time << "\t" << solutions[i].first.relMIPgap << "\t" << solutions[i].second.relMIPgap << endl;
		}

	}
	// close file
	out.close();
}