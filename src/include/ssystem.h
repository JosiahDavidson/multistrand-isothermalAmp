/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2024 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

/* SimulationSystem class header.
 * This is the main object which controls the entire simulated system. */

#ifndef __SSYSTEM_H__
#define __SSYSTEM_H__

#include <unordered_map>

#include "scomplexlist.h"
#include "simoptions.h"
#include "statespace.h"


typedef std::vector<bool> boolvector;
typedef std::vector<bool>::iterator boolvector_iterator;

namespace result_type {

const static std::string STR_ERROR = "error";
const static std::string STR_NAN = "nan";
const static std::string STR_NOINITIAL = "noinitial";
const static std::string STR_TIMEOUT = "timeout";

}

class SimulationSystem {
public:
	SimulationSystem(PyObject *options);
	SimulationSystem(PyObject* options, EnergyModel *energy_model);
	// helper method for constructors
	void parametrize(PyObject *ssystem);
	void parametrize(EnergyModel *energy_model);

	~SimulationSystem(void);

	void StartSimulation(void);
	void stateInfo(PyObject *start_state); // printing function
	void localTransitions(void); // builds all transitions in local statespace

	EnergyModel *GetEnergyModel(void);
	bool isFirstInstance(void);
	PyObject *calculateEnergy(PyObject *start_state, int typeflag);

private:
	void StartSimulation_Standard(void);
	void StartSimulation_FirstStep(void);
	void StartSimulation_Trajectory(void);
	void StartSimulation_Transition(void);

	void SimulationLoop_Standard(void);
	void SimulationLoop_FirstStep(void);
	void SimulationLoop_Trajectory(void);
	void SimulationLoop_Transition(void);

	int InitializeSystem(PyObject *alternate_start = NULL);

	void InitializePRNG();
	void generateNextRandom(void);
	void finalizeRun(void);
	void finalizeSimulation(void);

	// helper function for sending current state to Python side
	void dumpEndStateToPython(void);
	void sendTrajectory_CurrentStateToPython(double current_time, double arrType = -77.0);
	void sendTransitionStateVectorToPython(boolvector transition_states, double current_time);

	void exportTime(double& simTime, double& lastExportTime);
	void exportInterval(double simTime, long period, double arrType = -88.0);

	void initializeTransitions(PyObject*, bool);
	void iterateTransitions(
		PyObject *start_state, bool uni, bool send, bool print,
		std::vector<string>* uniMoveInfo = NULL);

	EnergyModel* energyModel = NULL;
	StrandComplex *startState = NULL;
	SComplexList *complexList = NULL;
	SimOptions *simOptions = NULL;

	PyObject *system_options = NULL;

	/* Management of the Libc PRNG buffer.
	 *
	 *   Initialization:
	 *     First trajectory:
	 *       - `SimulationSystem::InitializePRNG()`
	 *     Subsequent trajectories:
	 *       - `SimulationSystem::generateNextRandom()`
	 *
	 *   Read/write access during simulation:
	 *     Kinetic Monte Carlo:
	 *       - `SimTimer::advanceTime()`
	 *
	 *   Read access for simulation output:
	 *     Trajectory initialization:
	 *       - `SimulationSystem::InitializeSystem()`
	 *     Kinetic Monte Carlo:
	 *       - `SimulationSystem::sendTrajectory_CurrentStateToPython()`
	 *     Trajectory termination:
	 *       - `SimulationSystem::dumpEndStateToPython()`
	 *       - `SimOptions::stopResult*()`
	 */
	seed32_t trajectory_seed = 0; 		 // PRNG seed at the initial state
	seed48_t* current_state_seed = NULL; // full PRNG buffer at the current state

	long simulation_mode;
	long simulation_count_remaining;

	// triggers for output
	bool exportStatesTime = false;
	bool exportStatesInterval = false;

	// counters for timeouts and no-move initial states.
	int noInitialMoves = 0;
	int timeOut = 0;

	// A builder object that is only used if export is toggled
	Builder builder;

	bool firstInstance = false;

};

#endif
