/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2025 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

/* SimOption class header.
 * This is currently a wrapper for the options set in Python. */

#ifndef __SIMOPTIONS_H_
#define __SIMOPTIONS_H_


// Libc PRNG seed (32 bits) & buffer (48 bits) types
typedef long seed32_t;
typedef unsigned short seed48_t[3];


#include "energyoptions.h"

using std::vector;
using std::string;


struct complex_input {

	std::string sequence;
	std::string structure;
	identList* list;

	complex_input() { // empty constuctor
		sequence = "default";
		structure = "default";
		list = NULL;
	}

	complex_input(const char* string1, const char* string2, identList* list1) {

		char *tempseq = (char *) new char[strlen(string1) + 1];
		char *tempstruct = (char *) new char[strlen(string2) + 1];
		strcpy(tempseq, string1);
		strcpy(tempstruct, string2);

		sequence = string(tempseq);
		structure = string(tempstruct);
		list = list1;

	}
};


// FD: SimOptions contains an EnergyOptions object.
// EnergyModel contains all precomputed maps AND relevant energy functions.
//
// hiarchy: energy model > simoptions > energyoptions.

class SimOptions {
public:

	// Constructors
	virtual ~SimOptions(void) =0;

	// Non-virtual
	bool useTrajSeed();
	bool useStateSeed();
	long getTrajSeed();
	seed48_t* getStateSeed();
	EnergyOptions* getEnergyOptions();
	long getSimulationMode();
	long getSimulationCount();
	long getOInterval(void);
	double getOTime(void);
	long getStopOptions(void);
	long getStopCount(void);
	double getStartSimTime(void);
	double getMaxSimTime(void);

	bool getPrintIntialFirstStep(); // true if the initial state has to be exported.

	bool usingArrhenius(void);

	// Virtual methods

	virtual PyObject* getPythonSettings(void) = 0;
	virtual void generateComplexes(PyObject*, long) = 0;
	virtual stopComplexes* getStopComplexes(int) = 0;

	// Exit signalling
	virtual void stopResultError(long) = 0;
	virtual void stopResultNan(long) = 0;
	virtual void stopResultNormal(long, double, char*) = 0;
	virtual void stopResultTime(long, double) = 0;

	// For a first step result, also report the collision rate.
	virtual void stopResultFirstStep(long, double, double, const char*) = 0;


	// IO Methods
	string toString(void);

	// Cotranscriptional folding settings.
	bool cotranscriptional = false;
	double cotranscriptional_rate = 0.002; // delay between adding nucleotides (seconds)
	const int initialActiveNT = 8;	// initial number of active nucleotides.

	vector<complex_input>* myComplexes = NULL;
	EnergyOptions* energyOptions = NULL;

	bool statespaceActive;
	long verbosity = 1;
	double ms_version = 0.0;

	// switch for debugging output
	bool debug;

protected:

	long simulation_mode = 0;
	long simulation_count = 0;
	long o_interval = 0;
	double o_time = 0;
	long stop_options = 0;
	long stop_count = 0;
	double start_sim_time = 0;
	double max_sim_time = 0;

	seed32_t traj_seed = 0;
	seed48_t state_seed = {0, 0, 0};
	bool fixedTrajSeed = false;
	bool fixedStateSeed = false;

	stopComplexes* myStopComplexes = NULL;

	bool printInitialFirstStep = false;

};

class PSimOptions: public SimOptions {
public:
	//constructors
	PSimOptions(PyObject *options);

	static int checkPythonType(PyObject *options);
	static void clear(PyObject* options);
	PyObject* getPythonSettings(void);
	void generateComplexes(PyObject *alternate_start, seed32_t trajectory_seed);
	stopComplexes* getStopComplexes(int);

	// Error signaling
	void stopResultError(seed32_t);
	void stopResultNan(seed32_t);
	void stopResultNormal(seed32_t, double, char*);
	void stopResultTime(seed32_t, double);
	void stopResultFirstStep(seed32_t, double, double, const char*);

protected:
	PyObject *python_settings;

};

#endif
