/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2025 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#ifndef __SIMTIMER_H__
#define __SIMTIMER_H__

#include "simoptions.h"


class SimTimer {

public:
	SimTimer(SimOptions& myOptions);
	~SimTimer();

	void resetTime();
	void advanceTime();
	bool withinThresholds();
	bool wouldBeHit(const double);
	bool checkHit(const double);
	int checkHitBi(const double collisionRate);

	// initialize `SimTimer.seed`
	void setPRNG();
	// synchronize with Libc internal PRNG buffer
	static void readPRNG(seed48_t * prng);
	void writePRNG();
	// log PRNG state
	seed48_t* getPRNG();
	static void printPRNG(seed48_t * prng);

	bool checkForNewNucleotide(void);
	friend std::ostream& operator<<(std::ostream&, SimTimer&);

	double rate = 0.0;
	double stime = 0.0;
	long ssteps = 0;
	double startsimtime = 0.0;
	double maxsimtime = 0.0;
	long maxsimsteps = 0;
	double last_trajectory_time = 0.0;
	long stopcount = 0;
	long stopoptions = 0;

	int nuclAdded = 0;

	// inspection needs to set this to a non-random variable
	double rchoice = 0.0;

	// external buffer for Libc PRNG, used by `advanceTime()`
	seed48_t seed = {0, 0, 0};

	SimOptions* simOptions = NULL;

};

#endif
