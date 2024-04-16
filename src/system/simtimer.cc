/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2024 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#include "simtimer.h"


SimTimer::SimTimer(SimOptions& myOptions) {

	startsimtime = myOptions.getStartSimTime();
	maxsimtime = myOptions.getMaxSimTime();
	stopcount = myOptions.getStopCount();
	stopoptions = myOptions.getStopOptions();

	// saving the pointer to enable access to cotranscriptional timing values
	simOptions = &myOptions;

	stime = startsimtime;
	setPRNG();
}

// Populate `SimTimer.seed` from the Libc PRNG buffer, which was initialised
// earlier by `SimulationSystem.InitializePRNG()`.
void SimTimer::setPRNG() {

	readPRNG(&seed);

	// printPRNG(&seed);
}

// Read out the Libc internal PRNG buffer.
void SimTimer::readPRNG(seed48_t *prng) {

	static seed48_t zero_seed = {0, 0, 0};
	unsigned short* buf = seed48(zero_seed);
	memcpy(prng, buf, sizeof(*prng));
	seed48(*prng);
	buf = NULL;
}

// Overwrite the Libc internal PRNG buffer with `SimTimer.seed`. This is
// used to persist the PRNG state beyond the lifetime of this `SimTimer`
// instance, in particular for `SimulationSystem::generateNextRandom()`.
void SimTimer::writePRNG() {

	seed48(seed);
}

seed48_t* SimTimer::getPRNG() {

	return &seed;
}

void SimTimer::printPRNG(seed48_t *prng) {

	cout << "seed=("
		 << (*prng)[0] << "," << (*prng)[1] << "," << (*prng)[2]
		 << ")" << endl;
}

// Advance the simulation time according to the current `SimTimer.rate`.
// Conditioned on the random variables sampled here, the entire simulation
// is deterministic. Note that `SimTimer.seed` is used as an external
// seed buffer for the Libc PRNG.
void SimTimer::advanceTime() {

	rchoice = rate * erand48(seed);
	stime += log(1. / (1.0 - erand48(seed))) / rate;

	// printPRNG(&seed);
}

// returns TRUE if this transition needs to be executed.
bool SimTimer::wouldBeHit(const double rate) {

	return rchoice < rate;
}

// returns which Nth collision needs to be used
int SimTimer::checkHitBi(const double collisionRate) {

	return (int) floor(rchoice / collisionRate);

}

// returns TRUE if this transition needs to be executed
// and subtracts the rate
bool SimTimer::checkHit(const double rate) {

	rchoice = rchoice - rate;
	return (rchoice < 0);

}

// returns TRUE if a new nucleotide is to be added to the chain
bool SimTimer::checkForNewNucleotide(void) {

	if (simOptions->cotranscriptional && stime > (nuclAdded + simOptions->initialActiveNT) * simOptions->cotranscriptional_rate) {

		nuclAdded++;
		return true;
	}

	return false;
}

std::ostream& operator<<(std::ostream& ss, SimTimer& timer) {

	ss << "[rchoice=" << timer.rchoice
	   << ", rate=" << timer.rate
	   << ", simTime=" << timer.stime << "]" << std::endl;
	return ss;
}