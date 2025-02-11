/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2025 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#include <sstream>
#include <cstring>

#include "options.h"	 // Python options helper
#include "ssystem.h"

using std::vector;
using std::string;


SimOptions::~SimOptions(void) {}

PSimOptions::PSimOptions(PyObject* options) : python_settings(options)
{
	// initializers calling Python object -- these can use a super object getter.

	getBoolAttr(python_settings, initial_seed_flag, &fixedTrajSeed);
	if (fixedTrajSeed)
		getLongAttr(python_settings, initial_seed, &traj_seed);
	getBoolAttr(python_settings, state_seed_flag, &fixedStateSeed);
	if (fixedStateSeed)
		getPRNGSeed(python_settings, state_seed, state_seed);

	energyOptions = new PEnergyOptions(python_settings);

	getLongAttr(python_settings, simulation_mode, &simulation_mode);
	getLongAttr(python_settings, num_simulations, &simulation_count);
	getLongAttr(python_settings, output_interval, &o_interval);
	getDoubleAttr(python_settings, output_time, &o_time);
	getLongAttr(python_settings, stop_count, &stop_count);
	getLongAttr(python_settings, use_stop_conditions, &stop_options);
	getDoubleAttr(python_settings, simulation_start_time, &start_sim_time);
	getDoubleAttr(python_settings, simulation_time, &max_sim_time);

	// to print the initial state, used for statespace building
	getBoolAttr(python_settings, print_initial_first_step, &printInitialFirstStep);
	// to enable cotranscriptional folding, assumes a single strand
	getBoolAttr(python_settings, cotranscriptional, &cotranscriptional);
	getDoubleAttr(python_settings, cotranscriptional_rate, &cotranscriptional_rate);

	getLongAttr(python_settings, verbosity, &verbosity);
	getBoolAttr(python_settings, activestatespace, &statespaceActive);
	getDoubleAttr(python_settings, ms_version, &ms_version);

	// switch for debugging output
	debug = (verbosity > 2);

}

int PSimOptions::checkPythonType(PyObject* options) {
    if (strcmp(Py_TYPE(options)->tp_name, "Options") == 0)
        return 0;
    else {
        PyErr_SetString(PyExc_TypeError,
                        "Expected an argument of type `options.Options`.");
        return -1;
    }
}

void PSimOptions::clear(PyObject* options) {
	Py_DECREF(PyObject_CallMethod(options, "clear", "()"));
}

string SimOptions::toString() {

	std::stringstream ss;

	ss << "simulation_mode = " << simulation_mode << "\n";
	ss << "simulation_count = " << simulation_count << "\n";
	ss << "o_interval = " << o_interval << "\n";
	ss << "o_time = " << o_time << "\n";
	ss << "stop_options = " << stop_options << "\n";
	ss << "stop_count = " << stop_count << "\n";
	ss << "start_sim_time = " << start_sim_time << "\n";
	ss << "max_sim_time = " << max_sim_time << "\n";
	ss << "traj_seed = " << traj_seed << "\n";
	ss << "state_seed = ("
	   << state_seed[0] << "," << state_seed[1] << "," << state_seed[2]
	   << ")\n";

	string output = ss.str();
	output += energyOptions->toString();
	return output;

}

bool SimOptions::useTrajSeed() {

	return fixedTrajSeed;
}

bool SimOptions::useStateSeed() {

	return fixedStateSeed;
}

seed32_t SimOptions::getTrajSeed() {

	return traj_seed;
}

seed48_t* SimOptions::getStateSeed() {

	return &state_seed;
}

EnergyOptions* SimOptions::getEnergyOptions() {

	return energyOptions;

}

long SimOptions::getSimulationMode(void) {

	return simulation_mode;

}

long SimOptions::getSimulationCount(void) {

	return simulation_count;

}

long SimOptions::getOInterval(void) {

	return o_interval;

}

double SimOptions::getOTime(void) {

	return o_time;

}

long SimOptions::getStopOptions(void) {

	return stop_options;

}

long SimOptions::getStopCount(void) {

	return stop_count;

}

double SimOptions::getStartSimTime(void) {

	return start_sim_time;

}

double SimOptions::getMaxSimTime(void) {

	return max_sim_time;

}

bool SimOptions::getPrintIntialFirstStep() {

	return printInitialFirstStep;
}

bool SimOptions::usingArrhenius(void) {

	return energyOptions->usingArrhenius();

}

PyObject* PSimOptions::getPythonSettings() {

	return python_settings;

}

void PSimOptions::generateComplexes(PyObject *alternate_start, seed32_t trajectory_seed) {

	if (debug)
		cout << "Start generating complexes..." << endl;

	myComplexes = new vector<complex_input>(0); // wipe the pointer to the previous object;

	if (myComplexes->size() == 0) {
		const char *sequence, *structure;
		class identList *id;
		int start_count;
		PyObject *py_start_state = NULL, *py_complex = NULL;
		PyObject *py_seq = NULL, *py_struc = NULL;
		PyObject *py_err = NULL;

		if (alternate_start != NULL)
			py_start_state = alternate_start;
		else
			py_start_state = getListAttrReify(python_settings, start_state);

		if (py_start_state != Py_None)
			// doesn't need reference counting for this size call.
			// the getlistattr call we decref later.
			start_count = static_cast<int>(PyList_GET_SIZE(py_start_state));
		else
			start_count = 0;

		if (start_count == 0) {	// FD Jun 2018: adding throw if no initial state is set.
			throw std::invalid_argument("Initial state was not set.");
		}

		for (int index = 0; index < start_count; index++) {
			py_complex = PyList_GET_ITEM(py_start_state, index);
			// Borrowed reference, we do NOT decref it at end of loop.

#ifdef DEBUG_MACROS
			printPyError_withLineNumber();
#endif

			sequence = getStringAttrReify(py_complex, sequence, py_seq);
			// new reference

			structure = getStringAttrReify(py_complex, structure, py_struc);
			// new reference
			// Need to check if an error occurred, specifically, it could be an OSError due to sample failing. If so, we need to get the heck out of dodge right now.
			py_err = PyErr_Occurred();
			// py_err is a borrowed reference

			if (py_err != NULL) { // then an error occurred while getting the structure. Test for OSError (sample failure):
				if (PyErr_ExceptionMatches(PyExc_OSError)) {
					fprintf(stderr,
							"MULTISTRAND: Starting Structure could not be retrieved for index %d in your options object's start_state. This is likely due to Boltzmann sampling failing: please check that the program 'sample' exists and points correctly to the NUPACK sample binary. Or try 'print o.start_state[%d].structure' where 'o' is your options object and refer to that error message (if any).\n",
							index, index);
				} else {
					fprintf(stderr, "MULTISTRAND: An unidentified exception occurred while trying to initialize the system.\n");
				}
				return;
			}
			if (debug)
				printf("Complex %d: sequence='%s', structure='%s'\n", index, sequence, structure);

			id = getID_list(python_settings, index, alternate_start);

			complex_input myTempComplex = complex_input(sequence, structure, id);

			// StrandComplex does make its own copy of the seq/structure, so we can now decref.
			myComplexes->push_back(myTempComplex);

			Py_DECREF(py_seq);
			Py_DECREF(py_struc);
		}

		// Update the current seed and store the starting structures
		//   note: only if we actually have a system_options, e.g. no alternate start
		if (alternate_start == NULL && python_settings != NULL) {
			setLongAttr(python_settings, interface_trajectory_seed, trajectory_seed);
		}
		traj_seed = trajectory_seed;

	}

	return;
}

stopComplexes* PSimOptions::getStopComplexes(int) {

	myStopComplexes = getStopComplexList(python_settings, 0);

	return myStopComplexes;

}

void PSimOptions::stopResultError(seed32_t seed) {

	if (!statespaceActive) {
		printStatusLine(python_settings, seed, STOPRESULT_ERROR, 0.0, result_type::STR_ERROR.c_str());
	}

}

void PSimOptions::stopResultNan(seed32_t seed) {

	if (!statespaceActive) {
		printStatusLine(python_settings, seed, STOPRESULT_NAN, 0.0, result_type::STR_NAN.c_str());
	}

}

void PSimOptions::stopResultNormal(seed32_t seed, double time, char* message) {

	if (!statespaceActive) {
		printStatusLine(python_settings, seed, STOPRESULT_NORMAL, time, message);
	}

}

void PSimOptions::stopResultTime(seed32_t seed, double time) {

	if (!statespaceActive) {
		printStatusLine(python_settings, seed, STOPRESULT_TIME, time, result_type::STR_TIMEOUT.c_str());
	}

}

void PSimOptions::stopResultFirstStep(seed32_t seed, double stopTime, double rate, const char* message) {

	if (!statespaceActive) {
		printStatusLine_First_Bimolecular(python_settings, seed, STOPRESULT_NORMAL, stopTime, rate, message);
	}
}
