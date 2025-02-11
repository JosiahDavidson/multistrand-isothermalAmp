/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2025 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#include <time.h>

#include "options.h"
#include "multistrand_module.h"


SimulationSystem::SimulationSystem(PyObject *options) :
	system_options(options)
{
	simOptions = new PSimOptions(options);
    parametrize(SimSystemObject_lookup(options));
}

SimulationSystem::SimulationSystem(PyObject *options, EnergyModel *energy_model) :
	system_options(options)
{
	simOptions = new PSimOptions(options);
    parametrize(energy_model);
}

void SimulationSystem::parametrize(PyObject *ssystem) {
    // reuse `EnergyModel` if available
    EnergyModel *energy_model = NULL;
    if (ssystem != Py_None) {
        energy_model = ((SimSystemObject *) ssystem)->sys->GetEnergyModel();
    }
    parametrize(energy_model);
}

void SimulationSystem::parametrize(EnergyModel *energy_model) {
    simulation_mode = simOptions->getSimulationMode();
	simulation_count_remaining = simOptions->getSimulationCount();

    // reuse `EnergyModel` if available
    if (energy_model == NULL) {
		if (testLongAttr(system_options, parameter_type, =, 0))
			throw std::invalid_argument(
				"Attempting to load ViennaRNA parameters (deprecated)");
        energyModel = new NupackEnergyModel(simOptions->getPythonSettings());
        firstInstance = true;
	} else {
	    energyModel = energy_model;
	}

	exportStatesInterval = (simOptions->getOInterval() > 0);
	exportStatesTime = (simOptions->getOTime() >= 0);

	builder = Builder(simOptions);

	energyModel->writeConstantsToFile();

}

SimulationSystem::~SimulationSystem(void) {

	if (complexList != NULL) {
		delete complexList;
		complexList = NULL;
	}
	if (energyModel != NULL && firstInstance) {
		delete energyModel;
		energyModel = NULL;
	}

	// the remaining members are not our responsibility, we null them out
	// just in case something thread-unsafe happens.
	if (simOptions->myComplexes != NULL)
		delete simOptions->myComplexes;
	if (simOptions != NULL)
		delete simOptions;

}

void SimulationSystem::StartSimulation(void) {

	InitializePRNG();
	if (simulation_mode & SIMULATION_MODE_FLAG_FIRST_BIMOLECULAR) {
		StartSimulation_FirstStep();
	} else if (simulation_mode & SIMULATION_MODE_FLAG_TRAJECTORY) {
		StartSimulation_Trajectory();
	} else if (simulation_mode & SIMULATION_MODE_FLAG_TRANSITION) {
		StartSimulation_Transition();
	} else {
		StartSimulation_Standard();
	}
	finalizeSimulation();

}

EnergyModel *SimulationSystem::GetEnergyModel(void) {
    return energyModel;
}

bool SimulationSystem::isFirstInstance(void) {
    return firstInstance;
}

void SimulationSystem::StartSimulation_FirstStep(void) {

	while (simulation_count_remaining > 0) {
		if (InitializeSystem() != 0)
			return;

		SimulationLoop_FirstStep();
		finalizeRun();
	}

}

void SimulationSystem::StartSimulation_Standard(void) {

	while (simulation_count_remaining > 0) {
		if (InitializeSystem() != 0)
			return;

		SimulationLoop_Standard();
		finalizeRun();
	}

}

void SimulationSystem::StartSimulation_Transition(void) {

	while (simulation_count_remaining > 0) {
		if (InitializeSystem() != 0) {
			return;
		}

		SimulationLoop_Transition();
		finalizeRun();
	}
}

void SimulationSystem::StartSimulation_Trajectory(void) {

	while (simulation_count_remaining > 0) {

		if (InitializeSystem() != 0) {

			cout << "system not initialized; returning \n";
			return;
		}

		SimulationLoop_Trajectory();
		finalizeRun();
	}
}

void SimulationSystem::finalizeRun(void) {

	simulation_count_remaining--;
	pingAttr(system_options, increment_trajectory_count);

	generateNextRandom();

	// also ensure the builder does not remember the previous state
	builder.lastState = ExportData();

}

void SimulationSystem::finalizeSimulation(void) {

	if (simOptions->debug)
		cout << "Finalizing simulation..." << endl << flush;

	if (noInitialMoves > 0 and simOptions->verbosity) {

		cout << "No initial moves x" << noInitialMoves << "   (set verbosity = 0 to suppress) \n";
	}

	if (timeOut > 0 and simOptions->verbosity) {

		cout << "time-out detected x" << timeOut << "   (set verbosity = 0 to suppress) \n";

	}

	if (simOptions->statespaceActive) {

		builder.writeToFile();

	}

}

void SimulationSystem::SimulationLoop_Standard(void) {

	SimTimer myTimer(*simOptions);
	current_state_seed = myTimer.getPRNG();
	stopComplexes *traverse = NULL, *first = NULL;

	bool checkresult = false;

	complexList->initializeList();
	myTimer.rate = complexList->getTotalFlux();

	do {

		myTimer.advanceTime();

		if (myTimer.withinThresholds()) {
			// Why check here? Because we want to report the final state
			// as the one we were in before transitioning past the maximum
			// time, rather than the one after it. This is not entirely
			// obvious - to see why, look at a max time simulation for
			// finding an equilibrium distribution. What state are we
			// likely to observe after the time expires? Let's look at the
			// previous state - the more stable the state, the less total
			// rate out and thus the expected time for transition is
			// higher. So in our most likely case for what state we'll see
			// /after/ the max sim time is hit is actually the (more)
			// unstable state it transitioned to!

			// FD: Mathematically it is also the correct thing to do,
			// FD: when we remember the memoryless property of the Markov chain

			(void) complexList->doBasicChoice(myTimer);

			myTimer.rate = complexList->getTotalFlux();

			if (myTimer.stopoptions) {
				if (myTimer.stopcount <= 0) {
					simOptions->stopResultError(trajectory_seed);
					return;
				}

				checkresult = false;
				first = simOptions->getStopComplexes(0);
				checkresult = complexList->checkStopComplexList(first->citem);
				traverse = first;

				while (traverse->next != NULL && !checkresult) {
					traverse = traverse->next;
					checkresult = complexList->checkStopComplexList(traverse->citem);
				}
				// Note: we cannot delete first here if checkresult != 0,
				// as traverse->tag may be needed. It will get checked at
				// that point and deleted once traverse->tag is used,
				// later.
				if (!checkresult) {
					delete first;
				}
			}
		}
	} while (myTimer.withinThresholds() && !checkresult);

	if (myTimer.stime == NAN) {
		simOptions->stopResultNan(trajectory_seed);
	} else if (checkresult) {
		dumpEndStateToPython();
		simOptions->stopResultNormal(trajectory_seed, myTimer.stime, traverse->tag);
		delete first;
	} else { // stime >= maxsimtime ||  ssteps >= maxsimsteps
		dumpEndStateToPython();
		simOptions->stopResultTime(trajectory_seed, myTimer.maxsimtime);
	}

	current_state_seed = NULL;
}

void SimulationSystem::SimulationLoop_Trajectory() {

	SimTimer myTimer(*simOptions);
	current_state_seed = myTimer.getPRNG();
	stopComplexes *traverse = NULL, *first = NULL;

	bool stopFlag = false;
	long current_state_count = 0;

	complexList->initializeList();
	myTimer.rate = complexList->getTotalFlux();

	if (myTimer.stopoptions) {
		if (myTimer.stopcount <= 0) {
			simOptions->stopResultError(trajectory_seed);
			return;
		}
		first = simOptions->getStopComplexes(0);
	}

	// write the initial state:
	if (exportStatesInterval)
		exportInterval(myTimer.stime, current_state_count);

	do {

		myTimer.advanceTime();

		if (simOptions->debug) {
			cout << "Printing my complexlist! *************************************** \n";
			complexList->toString(cout);
			cout<< endl;
		}

		if (exportStatesTime) {
			exportTime(myTimer.stime, myTimer.last_trajectory_time);
		}

		double ArrMoveType = complexList->doBasicChoice(myTimer);
		myTimer.rate = complexList->getTotalFlux();
		current_state_count += 1;

		if (exportStatesInterval) {
			exportInterval(myTimer.stime, current_state_count, ArrMoveType);
		}

		if (myTimer.stopoptions) {

			stopFlag = false;
			stopFlag = complexList->checkStopComplexList(first->citem);
			traverse = first;
			while (traverse->next != NULL && !stopFlag) {
				traverse = traverse->next;
				stopFlag = complexList->checkStopComplexList(traverse->citem);
			}
		}

	} while (myTimer.withinThresholds() && !stopFlag);

	if (myTimer.stime == NAN) {
		simOptions->stopResultNan(trajectory_seed);
	} else if (stopFlag) {
		dumpEndStateToPython();
		simOptions->stopResultNormal(trajectory_seed, myTimer.stime, traverse->tag);
		// now export the tag to the builder as well
		builder.stopResultNormal(myTimer.stime, string(traverse->tag));
	} else {
		dumpEndStateToPython();
		simOptions->stopResultTime(trajectory_seed, myTimer.stime);
	}

	if (first != NULL)
		delete first;

	current_state_seed = NULL;
}

void SimulationSystem::SimulationLoop_Transition(void) {

	SimTimer myTimer(*simOptions);
	current_state_seed = myTimer.getPRNG();
	stopComplexes *traverse = NULL, *first = NULL;

	bool checkresult = false;
	bool stopFlag = false;
	bool state_changed = false;

	if (myTimer.stopcount <= 0 || !myTimer.stopoptions) {
		// this simulation mode MUST have some stop conditions set.
		simOptions->stopResultError(trajectory_seed);
		return;
	}

// figure out which stop entries should cause us to halt, update a bool vector to
// have true in the indices corresponding to which stop states are halting states.

	boolvector stop_entries;
	boolvector transition_states;
	stop_entries.resize(myTimer.stopcount, false);
	transition_states.resize(myTimer.stopcount, false);

	complexList->initializeList();

	first = simOptions->getStopComplexes(0);
	traverse = first;
	checkresult = false;

	for (int idx = 0; idx < myTimer.stopcount; idx++) {
		if (strstr(traverse->tag, "stop:") == traverse->tag)
			stop_entries[idx] = true;

		checkresult = complexList->checkStopComplexList(traverse->citem);

		transition_states[idx] = checkresult;
		traverse = traverse->next;
	}
	delete first;
	sendTransitionStateVectorToPython(transition_states, myTimer.stime);
// start

	myTimer.rate = complexList->getTotalFlux();
	state_changed = false;
	stopFlag = false;
	do {

		myTimer.advanceTime();

		if (myTimer.withinThresholds()) {
			// See note in SimulationLoop_Standard

			complexList->doBasicChoice(myTimer);
			myTimer.rate = complexList->getTotalFlux();

			// check if our transition state membership vector has changed
			first = simOptions->getStopComplexes(0);
			checkresult = false;
			traverse = first;

			for (int idx = 0; idx < myTimer.stopcount; idx++) {

				checkresult = complexList->checkStopComplexList(traverse->citem);

				if (checkresult && stop_entries[idx]) {
					// multiple stop states could suddenly be true, we add
					// a status line entry for the first one found.
					if (!stopFlag) {
						simOptions->stopResultNormal(trajectory_seed, myTimer.stime, traverse->tag);
					}

					stopFlag = true;
				}

				if (!state_changed && transition_states[idx] != checkresult) {
					state_changed = true;
				}

				transition_states[idx] = checkresult;
				traverse = traverse->next;
			}
			delete first; // we can do this now as we no longer need to
						  // save the stoplist until the loop exits, due
						  // to moving printStatusLine to immediately upon
						  // finding the stopping condition.
			if (state_changed) {
				sendTransitionStateVectorToPython(transition_states, myTimer.stime);
				state_changed = false;
			}
		}
	} while (myTimer.withinThresholds() && !stopFlag);

	if (myTimer.stime == NAN) {
		simOptions->stopResultNan(trajectory_seed);
	} else if (stopFlag) {
		dumpEndStateToPython();
		simOptions->stopResultNormal(trajectory_seed, myTimer.stime, traverse->tag);
	} else { // stime >= maxsimtime ||  ssteps >= maxsimsteps
		dumpEndStateToPython();
		simOptions->stopResultTime(trajectory_seed, myTimer.maxsimtime);
	}

	current_state_seed = NULL;
}

void SimulationSystem::SimulationLoop_FirstStep(void) {

	SimTimer myTimer(*simOptions);
	current_state_seed = myTimer.getPRNG();
	stopComplexes *traverse = NULL, *first = NULL;

	bool stopFlag = false;

	double frate = 0.0;

	long current_state_count = 0;

	complexList->initializeList();

	myTimer.rate = complexList->getJoinFlux();

	// if the toggle is set, export the initial state with arrType equal to flux
	// and timestamp -1
	if (simOptions->getPrintIntialFirstStep()) {
		exportInterval(-1.0, current_state_count, myTimer.rate);
	}

// scomplexlist returns a 0.0 rate if there was a single complex in
// the system, and a -1.0 rate if there are exactly 0 join moves. So
// the 0.0 rate should probably be caught, though if you use a
// single complex system for a starting state it's probably
// deserved.

	if (myTimer.rate == 0.0) { // no initial moves

		noInitialMoves++;

		simOptions->stopResultFirstStep(trajectory_seed, 0.0, 0.0, result_type::STR_NOINITIAL.c_str());
		return;
	}

	myTimer.advanceTime(); // select an rchoice
	myTimer.startsimtime = 0.0; // but reset the jump in time, because first step mode.
	myTimer.resetTime();

	double ArrMoveType = complexList->doJoinChoice(myTimer);

	if (exportStatesInterval) {
		exportInterval(myTimer.stime, current_state_count, ArrMoveType);
	}

// store the forward rate used for the initial step so we can record it.
	frate = myTimer.rate * energyModel->getJoinRate_NoVolumeTerm() / energyModel->getJoinRate();

// This join rate is the dG_volume * bimolecular scaling constant
// used for forward transitions.  What we actually need is the
// bimolecular scaling constant * total move count. (dG volume is
// the volume dependent term that is not actually related to the
// 'collision' rate, but rather the volume we are simulating.

// Begin normal steps.
	myTimer.rate = complexList->getTotalFlux();

	do {

		myTimer.advanceTime();

		if (simOptions->debug) {
			cout << "Printing my complexlist! *************************************** \n";
			complexList->toString(cout);
			cout << endl;
		}

		// trajectory output via outputtime option
		if (exportStatesTime) {
			exportTime(myTimer.stime, myTimer.last_trajectory_time);
		}

		double ArrMoveType = complexList->doBasicChoice(myTimer);

		myTimer.rate = complexList->getTotalFlux();
		current_state_count++;

		if (exportStatesInterval) {
			exportInterval(myTimer.stime, current_state_count, ArrMoveType);
		}

		if (myTimer.stopcount > 0 && myTimer.stopoptions) {

			stopFlag = false;
			first = simOptions->getStopComplexes(0);
			traverse = first;
			stopFlag = complexList->checkStopComplexList(traverse->citem);

			while (traverse->next != NULL && !stopFlag) {
				traverse = traverse->next;
				stopFlag = complexList->checkStopComplexList(traverse->citem);
			}

			if (!stopFlag && first != NULL)
				delete first;
		}

	} while (myTimer.withinThresholds() && !stopFlag);

	if (stopFlag) {
		dumpEndStateToPython();
		simOptions->stopResultFirstStep(trajectory_seed, myTimer.stime, frate, traverse->tag);
		delete first;
	} else {
		timeOut++;
		dumpEndStateToPython();
		simOptions->stopResultFirstStep(trajectory_seed, myTimer.stime, frate, result_type::STR_TIMEOUT.c_str());
	}

	current_state_seed = NULL;
}

// Helper function to send the current state to the Python module.
void SimulationSystem::dumpEndStateToPython(void) {

	SComplexListEntry *temp = complexList->getFirst();
	ExportData data;

	while (temp != NULL) {
		temp->dumpComplexEntryToPython(data, energyModel);
		printComplexStateLine(simOptions->getPythonSettings(), current_state_seed, data);
		temp = temp->next;
	}
}

// Helper function to prepare a Python list object containing the bool information
// about which transition states we are in.
void SimulationSystem::sendTransitionStateVectorToPython(boolvector transition_states, double current_time) {

	PyObject *mylist = PyList_New((Py_ssize_t) transition_states.size());
	// we now have a new reference here that we'll need to DECREF.

	if (mylist == NULL)
		return; // TODO: Perhaps raise an exception to the Python side here that
				// we couldn't pass the information back...

	boolvector_iterator it;
	Py_ssize_t index = 0;

	for (it = transition_states.begin(); it < transition_states.end(); it++) {
		if (*it) // bool value was true
		{
			Py_INCREF(Py_True);
			PyList_SET_ITEM(mylist, index, Py_True);
			// ownership of this reference to Py_True has now been stolen by PyList_SET_ITEM.
		} else {
			Py_INCREF(Py_False);
			PyList_SET_ITEM(mylist, index, Py_False);
			// ownership of this reference to Py_False has now been stolen by PyList_SET_ITEM.
		}
		index++;
	}

	PyObject *transition_tuple = Py_BuildValue("dO", current_time, mylist);
	// we now have a new reference to transition_tuple.  note that our
	// reference to mylist has NOT been stolen by this call [it was
	// increffed in the call itself, so we can decref now with no
	// worries]
	Py_DECREF(mylist);
	// we now have no references to mylist directly owned [though transition tuple has one via BuildValue]

	pushTransitionInfo(system_options, transition_tuple);

	// transition_tuple has been decreffed by this macro, so we no longer own any references to it
}

// Helper function to send the current state to the Python module.
void SimulationSystem::sendTrajectory_CurrentStateToPython(double current_time, double arrType) {

	ExportData data;
	ExportData mergedData;

	SComplexListEntry *temp = complexList->getFirst();

	while (temp != NULL) {

		temp->dumpComplexEntryToPython(data, energyModel);
		if (!simOptions->statespaceActive)
			pushTrajectoryComplex(system_options, current_state_seed, data);
		mergedData.merge(data);
		temp = temp->next;
	}

	// for now, keep exporting the state to the regular interface too.
	if (simOptions->statespaceActive) {
		builder.addState(mergedData, arrType);
	} else {
		pushTrajectoryInfo(system_options, current_time);
		pushTrajectoryInfo2(system_options, arrType);
	}

}

// FD: OK to have alternate_start = NULL
int SimulationSystem::InitializeSystem(PyObject *alternate_start) {
	StrandComplex *tempcomplex;
	identList *id;

	simOptions->generateComplexes(alternate_start, trajectory_seed);

	// FD: Somehow, check if complex list is pre-populated.
	startState = NULL;
	if (complexList != NULL)
		delete complexList;

	if (simOptions->debug)
		cout << "myComplexes.size = " << simOptions->myComplexes->size() << endl << flush;

	complexList = new SComplexList(energyModel);

	// FD: this is the python - C interface
	for (unsigned int i = 0; i < simOptions->myComplexes->size(); i++) {

		char* tempSequence = utility::copyToCharArray(simOptions->myComplexes->at(i).sequence);
		char* tempStructure = utility::copyToCharArray(simOptions->myComplexes->at(i).structure);

		id = simOptions->myComplexes->at(i).list;

		tempcomplex = new StrandComplex(tempSequence, tempStructure, id, simOptions->debug);

		startState = tempcomplex;
		complexList->addComplex(tempcomplex);
	}

	if (simOptions->debug)
		cout << "Done initializing!" << endl;

	return 0;
}

void SimulationSystem::InitializePRNG() {

	if (simOptions->useStateSeed()) {

		trajectory_seed = 0;
		current_state_seed = simOptions->getStateSeed();
		// initialize the Libc internal PRNG buffer
		seed48(*current_state_seed);

	} else if (simOptions->useTrajSeed()) {

		trajectory_seed = simOptions->getTrajSeed();
		// initialize the Libc internal PRNG buffer
		srand48(trajectory_seed);

	} else {

		FILE *fp = NULL;
		if ((fp = fopen("/dev/urandom", "r")) != NULL) {
			// if urandom exists, use it to provide a seed
			long deviceseed;
			(void) fread(&deviceseed, sizeof(long), 1, fp);
			trajectory_seed = deviceseed;
			fclose(fp);
		} else {
			// use the possibly flawed time as a seed
			trajectory_seed = time(NULL);
		}
		// initialize the Libc internal PRNG buffer
		srand48(trajectory_seed);

	}

	current_state_seed = NULL;
}

void SimulationSystem::generateNextRandom(void) {

	seed48_t tmp_seed;
	if (simOptions->debug) {
		cout << "generateNextRandom:" << endl;
		SimTimer::readPRNG(&tmp_seed);
		cout << "  before: "; SimTimer::printPRNG(&tmp_seed);
	}

	// Generate the initial seed for the next trajectory from the last
	// PRNG state of the current trajectory. Assumes that `myTimer.writePRNG()`
	// was called in `myTimer.~SimTimer()` at the end of the current trajectory.
	trajectory_seed = lrand48();

	// initialize the Libc internal PRNG buffer
	srand48(trajectory_seed);
	current_state_seed = NULL;

	if (simOptions->debug) {
		SimTimer::readPRNG(&tmp_seed);
		cout << "  after: "; SimTimer::printPRNG(&tmp_seed);
	}
}

PyObject *SimulationSystem::calculateEnergy(PyObject *start_state, int typeflag) {
	double *values = NULL;
	PyObject *retval = NULL;

// calc based on current state, do not clean up anything.
	if (start_state != Py_None) {
		InitializeSystem(start_state);
		complexList->initializeList();
	}

	values = complexList->getEnergy(typeflag); // NUPACK energy output : bimolecular penalty, no Volume term.

	if (simOptions->debug)
		cout << "Done calculating energy!" << endl << flush;

	retval = PyTuple_New(complexList->getCount());
// New Reference, we return it.
// The complex list is a linked list and new items are added at the head; so we need to reverse the resulting list to get the data back out.
	for (int loop = complexList->getCount() - 1; loop >= 0; loop--)
		PyTuple_SET_ITEM(retval, loop, PyFloat_FromDouble(values[loop]));
// the reference from PyFloat_FromDouble is immediately stolen by PyTuple_SET_ITEM.

	delete[] values;
	return retval;
}

void SimulationSystem::exportTime(double& simTime, double& lastExportTime) {

	if (simTime - lastExportTime > simOptions->getOTime()) {
		lastExportTime += simOptions->getOTime();
		sendTrajectory_CurrentStateToPython(lastExportTime);
	}

}

void SimulationSystem::exportInterval(double simTime, long transitionCount, double arrType) {

	if ((transitionCount % simOptions->getOInterval()) == 0)
		sendTrajectory_CurrentStateToPython(simTime, arrType);

}

/*
 * uniform=false: correct elementary rates
 * uniform=true: uniform event space parametrization (all elementary rates = 1.0)
 */
void SimulationSystem::initializeTransitions(PyObject *start_state, bool uniform) {

	InitializeSystem(start_state);
	energyModel->inspection = uniform;

	complexList->initializeList();
	// required to compute join rates
	complexList->updateOpenInfo();
	// required to set cumulative join rate
	complexList->getTotalFlux();

}

/*
 * Given string representations of sequence and structure (held by a Python object),
 * wastefully iterate through adjacent states, i.e., obtain the (C++) `StrandComplex`
 * instances resulting from applying `StrandComplex::doChoice()` for each `Move`
 * in the `MoveContainer`. In mathematical terms, `StrandComplex` is a state representation,
 * whereas `Move` is a state action representation.
 *
 * This recursion scheme was introduced for `SimulationSystem::localTransitions()`
 * in the context of the "pathway elaboration" method, and is now also used in refactored
 * form to add adjacent state information to the state introspection in
 * `SimulationSystem::stateInfo()`. These two evaluation strategies are controlled
 * by the `send` and `print` arguments.
 * 
 * NOTE:
 * -----
 *
 * Unfortunately, extracting this information currently requires reinitializing the
 * `SimulationSystem` *for each move*, i.e., repeatedly recomputing the entire loop graph
 * and the associated tree of cumulative transition rates, even though the initial state
 * is semantically constant. The state action is then executed for each move by hijacking
 * the kinetic Monte Carlo (SSA) machinery as a random variable over elementary transitions,
 * and the outer loop uses `SimTimer` as a parametrization of this random variable's
 * event space, after artificially uniformizing all elementary rates.
 *
 * This workaround is a consequence of Multistrand's core design, in which the runtime
 * representation of states and transitions is a network of mutable objects with mutual
 * pointer references, and in which state actions are implemented as a wave of object
 * creations, mutations and deletions in this graph. Hence, the information required
 * to recover a previous state is in part deleted and in part spread over the heap, and
 * `InitializeSystem()` ends up being the least redundant entry point in the existing
 * code base for safely reconstructing a previous state.
 
 * If state actions (and their pre- and post-conditions) had *serializable specifications*
 * instead, -- as is already the case for states via the dot-parentheses notation --
 * then the simulation loop could be implemented as a "handler" for "transactions" or
 * "algebraic effects", and this method could use reverse actions (or deep copies)
 * to achieve its task more efficiently. This, however, would amount to a reimplementation
 * of all core data structures.
 */
void SimulationSystem::iterateTransitions(
	PyObject *start_state, bool uni, bool send, bool print,
	std::vector<string>* uniMoveInfo) {

	initializeTransitions(start_state, true);

	// iterator bounds
	int Nuni = complexList->getMoveCount(),
	    Nbi = int(round(complexList->getJoinFlux()));
	int lo = uni ? Nbi : 0,
	    hi = uni ? Nbi + Nuni : Nbi;

	string csep, lsep;
	if (print) {
	    csep  = "----------------------------------------";
	    lsep  = "********************";
		cout << (uni ? "[uni]" : "[bi]") << endl;
	}

	// iterator state
	int complexId0 = 0, complexId1 = 0, complexId0_new = 0, complexId1_new = 0,
	    loopIdx = 0, loopIdx_new = 0;
	bool complexChange;
	std::vector<string>::iterator uniMoveIter;
	if (uniMoveInfo != NULL) {
		assert (print);
		assert (uni);
		uniMoveIter = uniMoveInfo->begin();
	}
	stringstream info0, info1;

	// iteration
	for (int i = lo; i < hi; i++) {
		
		if (i > lo)
			initializeTransitions(start_state, true);

		// rely on event space parametrization
		SimTimer myTimer(*simOptions);
		current_state_seed = myTimer.getPRNG();
		myTimer.rchoice = i + 0.01;

		// export the initial state
		if (send && exportStatesInterval)
			exportInterval(myTimer.stime, 0, 8888);

		// move the state using the i-th transition
		double ArrMoveType = complexList->doBasicChoice(
			myTimer, &complexId0_new, &complexId1_new, &loopIdx_new,
			&info0, &info1);

		if (print) {
			// print complex/complex pair/loop subheaders on change
			complexChange =
				i == lo || complexId0_new != complexId0
				        || complexId1_new != complexId1;
			if (complexChange) {
				cout << csep << endl;
				if (uni)
					cout << "Complex" << complexId0_new << endl;
				else
					cout << "Complexes: "
					     << complexId0_new << "," << complexId1_new << endl
					     << info0.str();
				cout << csep << endl;
			}
			if (complexChange || loopIdx_new != loopIdx) {
				if (uni)
					cout << lsep << endl
						 << info0.str()
						 << lsep << endl;
			}

			// print move/join information
			if (uni)
				cout << *uniMoveIter;
			else
				cout << "Join" << i << " | " << info1.str() << endl;

			// print adjacent state information
			complexList->printStructures(cout);

			// advance iterator
			complexId0 = complexId0_new;
			complexId1 = complexId1_new;
			complexId0_new = complexId1_new = 0;
			loopIdx = loopIdx_new;
			loopIdx_new = 0;
			if (uni)
				uniMoveIter++;
			else
				stringstream().swap(info1);
			stringstream().swap(info0);
		}

		// export the state after the transition
		if (send) {
			if (exportStatesInterval)
				exportInterval(myTimer.stime, 1, ArrMoveType);
			finalizeRun();
		}
		current_state_seed = NULL;
	}
	
	if (print)
		cout << endl << flush;
}

void SimulationSystem::stateInfo(PyObject *start_state) {

	initializeTransitions(start_state, false);

	// print strand complexes
	string hline = "========================================";
	cout << hline << endl;
	complexList->printComplexList(cout);

	// print (cumulative) transition rates
	int moveCount = complexList->getMoveCount();
	double moveFlux = complexList->getMoveFlux(),
	       joinFlux = complexList->getJoinFlux(),
		   joinRate = energyModel->applyPrefactors(
	           energyModel->getJoinRate(), loopMove, loopMove),
		   joinConc = energyModel->simOptions->energyOptions->getJoinConcentration();
	cout << endl << hline << endl
	     << "# unimol. steps   = " << moveCount << endl
		 // this ratio assumes a Metropolis model (context-independent rates)
	     << "# bimol.  steps   = " << int(joinFlux / joinRate) << endl
		 << "unimol. cum. rate = " << moveFlux << " /s" << endl
		 << "bimol.  cum. rate = " << joinFlux << " /s" << endl
	     << "bimol.  step rate = " << joinRate << " /s" << endl
	     << "join conc.        = " << joinConc << endl
		 << endl << flush;
	
	// store unimolecular elementary transitions
	// (correct rates, but without adjacent states)
	std::vector<string> uniMoveInfo;
	uniMoveInfo.reserve(moveCount);
	complexList->printAllMoves(uniMoveInfo);

	// print elementary moves & adjacent states
	// (iterate with incorrect rates, splice in correct rates)
	bool send = false, print = true;
	cout << hline << endl;
	iterateTransitions(start_state, true, send, print, &uniMoveInfo);
	cout << hline << endl;
	iterateTransitions(start_state, false, send, print);

}

// FD: Build the statespace by taking N transitions, where N is the number of transitions out of the state
// If you call this function, the system will assume the active statespace toggle is set.
void SimulationSystem::localTransitions(void) {

	assert(simOptions->statespaceActive);
	InitializePRNG(); // the output dir will be '0' if unset

	bool send = true, print = false;
	iterateTransitions(NULL, false, send, print);
	iterateTransitions(NULL, true, send, print);

	finalizeSimulation();

}
