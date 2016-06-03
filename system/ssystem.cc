/*
  Copyright (c) 2007-2010 Caltech. All rights reserved.
  Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/

#include "options.h"
#include "ssystem.h"

#include <string.h>
#include <time.h>
#include <stdlib.h>

#ifdef PROFILING
#include "google/profiler.h"
#include "google/heap-profiler.h"
#endif

SimulationSystem::SimulationSystem( int argc, char **argv )
{
  return;
  // entire function is FAIL. Needs replacing.
}



SimulationSystem::SimulationSystem( PyObject *system_o )
{
  bool hflag = false;
#ifdef PROFILING
  if (!IsHeapProfilerRunning())
    {
      HeapProfilerStart("ssystem_init.heap");
      hflag = true;
    }
  ProfilerStart("ssystem_init_profile.prof");
#endif

  system_options = system_o;
  // We no longer need the below line; we are guaranteed that options
  // will have a good reference for the lifetime of our object, as the
  // controlling wrapper in multistrand_module.cc grabs the reference.

  //Py_INCREF( system_options );

  getLongAttr(system_options, simulation_mode, &simulation_mode );
  getLongAttr(system_options, num_simulations, &simulation_count_remaining);
  if( Loop::GetEnergyModel() == NULL)
    {
      dnaEnergyModel = NULL;

      if(  testLongAttr(system_options, parameter_type,=,0) )
        dnaEnergyModel = new ViennaEnergyModel( system_options );
      else
        dnaEnergyModel = new NupackEnergyModel( system_options );
      Loop::SetEnergyModel( dnaEnergyModel );
    }
  else
    {
      dnaEnergyModel = Loop::GetEnergyModel();
    }
  
  startState = NULL;
  complexList = NULL; 
#ifdef PROFILING
  ProfilerStop();
  if (hflag)
    {
    HeapProfilerDump("init");
    HeapProfilerStop();
    }
#endif
}

SimulationSystem::SimulationSystem( void )
{
  simulation_mode = -1;
  simulation_count_remaining = -1;

  if( Loop::GetEnergyModel() == NULL)
    {
      dnaEnergyModel = NULL;
    }
  else
    {
      dnaEnergyModel = Loop::GetEnergyModel();
    }

  system_options = NULL;
  startState = NULL;
  complexList = NULL; 
}

int SimulationSystem::getErrorFlag( void )
{
  if( dnaEnergyModel == NULL )
    return 1;
  return 0;
}

SimulationSystem::~SimulationSystem( void )
{
  if( complexList != NULL )
    delete complexList;  
  complexList = NULL;

  // the remaining members are not our responsibility, we null them out
  // just in case something thread-unsafe happens.
  dnaEnergyModel = NULL;

  // now handled in multistrand_module.cc:
  //  if( system_options != NULL )
  //  Py_DECREF( system_options );

  system_options = NULL;
  startState = NULL;
}


void SimulationSystem::InfoInitial( void )
{
  bool hflag = false;

  printf("Starting information dump about the initial state.	\n");
  printf("Initializing system now \n");

  if( InitializeSystem() != 0)
	return;

  printf("Initializing list now \n");
  complexList->initializeList();

  double rate =  complexList->getJoinFlux();
  printf("The joinFlux is %f \n", rate);

  rate = complexList->getTotalFlux();
  printf("The totalFlux is %f \n", rate);



}


void SimulationSystem::StartSimulation( void )
{
  bool hflag = false;
#ifdef PROFILING
  if (!IsHeapProfilerRunning())
    {
      HeapProfilerStart("ssystem_run.heap");
      hflag = true;
    }
  ProfilerStart("ssystem_run_profile.prof");
#endif
  printf("Simulation Mode TEST: %d\n",simulation_mode);
  if( simulation_mode & SIMULATION_MODE_FLAG_FIRST_BIMOLECULAR )
    {
      StartSimulation_FirstStep();
    }
  else if( simulation_mode & SIMULATION_MODE_FLAG_TRAJECTORY )
    {
      StartSimulation_Trajectory();
    }
  else if( simulation_mode & SIMULATION_MODE_FLAG_TRANSITION )
    {
      StartSimulation_Transition();
    }
  else
    StartSimulation_Standard();

#ifdef PROFILING
  ProfilerStop();
  if (hflag)
    {
      HeapProfilerDump("final");
      HeapProfilerStop();
    }
#endif
}


void SimulationSystem::StartSimulation_Standard( void )
{
  InitializeRNG();  
  while( simulation_count_remaining > 0)
    {
      if( InitializeSystem() != 0)
		return;

      SimulationLoop_Standard();

      simulation_count_remaining--;
      pingAttr( system_options, increment_trajectory_count );

      generateNextRandom();
    }
}

void SimulationSystem::StartSimulation_Transition( void )
{
  InitializeRNG();  
  while( simulation_count_remaining > 0)
    {
      if( InitializeSystem() != 0)
		return;

      SimulationLoop_Transition();

      simulation_count_remaining--;
      pingAttr( system_options, increment_trajectory_count );
      generateNextRandom();
    }
}


void SimulationSystem::StartSimulation_Trajectory( void )
{
  long ointerval;
  double otime;

  getLongAttr(system_options, output_interval,&ointerval);
  getDoubleAttr(system_options, output_time,&otime);

  InitializeRNG();  
  while( simulation_count_remaining > 0)
    {
      if( InitializeSystem() != 0)
		return;

      SimulationLoop_Trajectory( ointerval, otime );

      simulation_count_remaining--;
      pingAttr( system_options, increment_trajectory_count );

      generateNextRandom();
    }
}


void SimulationSystem::SimulationLoop_Standard( void )
{
  double rchoice,rate,stime,ctime;
  // Could really use some commenting on these local vars.
  rchoice = rate = stime = ctime = 0.0;


  double maxsimtime;
  maxsimtime = -1.0;
  
  int curcount = 0;
  bool checkresult = false;
  long stopcount = 0, stopoptions = 0;
  class stopcomplexes *traverse = NULL, *first=NULL;

  getLongAttr(system_options, use_stop_conditions,&stopoptions);
  getLongAttr(system_options, stop_count,&stopcount);
  getDoubleAttr(system_options, simulation_time,&maxsimtime);

  complexList->initializeList();

  rate = complexList->getTotalFlux();
 
  do {

    rchoice = rate * drand48();
    
    stime += (log( 1. / (1.0 - drand48()) ) / rate ); 
    // 1.0 - drand as drand returns in the [0.0, 1.0) range, we need a (0.0,1.0] range.
    // see notes below in First Step mode.

	if( stime < maxsimtime )  
	  {
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

		// FD: Mathematically, this is also the right move, because
		// the system is in the previous state when the time-bound passed.
		// See also the memoryless property of Markov chains.

		complexList->doBasicChoice( rchoice, stime );
		rate = complexList->getTotalFlux();

		if( stopoptions )
		  {
			if( stopcount <= 0 )
			  {
				printStatusLine(system_options, current_seed, STOPRESULT_ERROR, 0.0, NULL);
				return; 
			  }
			checkresult = false;
			first = getStopComplexList( system_options, 0 );
			checkresult = complexList->checkStopComplexList( first->citem );
			traverse = first;
			while( traverse->next != NULL && !checkresult )
			  {
				traverse = traverse->next;
				checkresult = complexList->checkStopComplexList( traverse->citem );
			  }
			// Note: we cannot delete first here if checkresult != 0,
			// as traverse->tag may be needed. It will get checked at
			// that point and deleted once traverse->tag is used,
			// later.
			if( !checkresult )
			  delete first;
		  }
	  }
  } while( stime < maxsimtime && !checkresult);
  
  if( stime == NAN )
    printStatusLine(system_options, current_seed, STOPRESULT_NAN, 0.0, NULL);
  else if ( checkresult  )
    {
	  dumpCurrentStateToPython();
      printStatusLine(system_options, current_seed, STOPRESULT_NORMAL, stime, traverse->tag );
      delete first;
    }
  else // stime >= maxsimtime
	{
	  dumpCurrentStateToPython();
	  printStatusLine(system_options, current_seed, STOPRESULT_TIME, maxsimtime, NULL );
	}
}

void SimulationSystem::SimulationLoop_Trajectory( long output_count_interval, double output_time_interval )
{
  double rchoice, rate, current_simulation_time, last_trajectory_time, maxsimtime;
  // Could really use some commenting on these local vars.
  rchoice = rate = 0.0;
  maxsimtime = -1.0;

  bool checkresult = false;
  long current_state_count = 0;
  long stopcount = 0, stopoptions = 0;
  class stopcomplexes *traverse = NULL, *first=NULL;

  // The tag of the stop state reached. 
  //char *tag = NULL;
  
  //  SComplexListEntry *affectedComplex = NULL;

  getDoubleAttr(system_options, simulation_time,&maxsimtime);

  complexList->initializeList();
  rate = complexList->getTotalFlux();
  
  // We start at the beginning of time.
  current_simulation_time = 0.0;
  
  // The last time we gave the output state.
  last_trajectory_time = 0.0;


  getLongAttr(system_options, use_stop_conditions,&stopoptions);
  getLongAttr(system_options, stop_count,&stopcount);

  if( stopoptions )
	{
	  if( stopcount <= 0 )
		{
		  printStatusLine(system_options, current_seed, STOPRESULT_ERROR, 0.0, NULL);
		  return; 
		}
	  first = getStopComplexList( system_options, 0 );
	}


  do {
    rchoice = rate * drand48();

    current_simulation_time += (log( 1. / (1.0 - drand48()) ) / rate ); 

    // trajectory output via outputtime option
	// we check this here so the reported state is the one present at the time
	// listed, rather than the one /after/ that.
    if( output_time_interval > 0.0 )
      if( current_simulation_time - last_trajectory_time > output_time_interval )
        {
          last_trajectory_time += output_time_interval;
		  sendTrajectory_CurrentStateToPython( last_trajectory_time );
          //complexList->printComplexList( 0 );
        }


    // 1.0 - drand as drand returns in the [0.0, 1.0) range, we need a (0.0,1.0] range.
    // see notes below in First Step mode.

    complexList->doBasicChoice( rchoice, current_simulation_time );
    rate = complexList->getTotalFlux();
    current_state_count += 1;
	

	if( stopoptions )
	  {
		checkresult = false;
		checkresult = complexList->checkStopComplexList( first->citem );
		traverse = first;
		while( traverse->next != NULL && !checkresult )
		  {
			traverse = traverse->next;
			checkresult = complexList->checkStopComplexList( traverse->citem );
		  }
	  }

    //    trajectory output via outputinterval option
    if( output_count_interval >= 0 )
      if( (current_state_count % output_count_interval) == 0)
		{
		  sendTrajectory_CurrentStateToPython( current_simulation_time );
		  //         complexList->printComplexList( 0 );
		}


  } while( current_simulation_time < maxsimtime && !checkresult );
            
  if( current_simulation_time == NAN )
    printStatusLine(system_options, current_seed, 
                    STOPRESULT_NAN, 0.0, 
                    NULL);
  else if ( checkresult )
    {
      printStatusLine(system_options,    current_seed, 
                      STOPRESULT_NORMAL, current_simulation_time, 
                      traverse->tag );
    }
  else
    printStatusLine(system_options,  current_seed, 
                    STOPRESULT_TIME, current_simulation_time, 
                    NULL );

  if( first != NULL )
	delete first;
}

void SimulationSystem::SimulationLoop_Transition( void )
{
  double rchoice,rate,stime,ctime;
  // Could really use some commenting on these local vars.
  rchoice = rate = stime = ctime = 0.0;

  double maxsimtime, otime;
  maxsimtime = otime = -1.0;
  
  bool checkresult = false;
  bool stop_flag = false;
  bool state_changed = false;
  long stopcount = 0, stopoptions = 0, sMode = 0;
  long ointerval = -1;
  class stopcomplexes *traverse = NULL, *first=NULL;

  getLongAttr(system_options, simulation_mode,&sMode); 
  getLongAttr(system_options, output_interval,&ointerval);
  getLongAttr(system_options, use_stop_conditions,&stopoptions);
  getLongAttr(system_options, stop_count,&stopcount);
  getDoubleAttr(system_options, simulation_time,&maxsimtime);
  getDoubleAttr(system_options, output_time,&otime);

  if( stopcount <= 0 || !stopoptions)
	{
	  // this simulation mode MUST have some stop conditions set.
	  printStatusLine(system_options, current_seed, STOPRESULT_ERROR, 0.0, NULL);
	  return; 
	}


  // figure out which stop entries should cause us to halt, update a bool vector to
  // have true in the indices corresponding to which stop states are halting states.

  boolvector stop_entries;
  boolvector transition_states;
  stop_entries.resize( stopcount, false );
  transition_states.resize( stopcount, false );

  complexList->initializeList();

  first = getStopComplexList( system_options, 0 );
  traverse = first;
  checkresult = false;
  for( int idx = 0; idx < stopcount; idx++ )
	{
	  if( strstr( traverse->tag, "stop:") == traverse->tag )
		stop_entries[idx] = true;

	  checkresult = complexList->checkStopComplexList( traverse->citem );

	  transition_states[idx] = checkresult;
	  traverse = traverse->next;
	}
  delete first;
  sendTransitionStateVectorToPython( transition_states, stime );
  // start


  rate = complexList->getTotalFlux();
  state_changed = false;
  stop_flag = false;
  do {

    rchoice = rate * drand48();
    
    stime += (log( 1. / (1.0 - drand48()) ) / rate ); 
    // 1.0 - drand as drand returns in the [0.0, 1.0) range, we need a (0.0,1.0] range.
    // see notes below in First Step mode.

	if( stime < maxsimtime )  
	  {
		// See note in SimulationLoop_Standard

		complexList->doBasicChoice( rchoice, stime );
		rate = complexList->getTotalFlux();

		// check if our transition state membership vector has changed
		checkresult = false;
		first = getStopComplexList( system_options, 0 );
		traverse = first;
		for( int idx = 0; idx < stopcount; idx++ )
		  {
			checkresult = complexList->checkStopComplexList( traverse->citem );
			if( checkresult && stop_entries[idx] == true)
			  {
				// multiple stop states could suddenly be true, we add
				// a status line entry for the first one found.
				if( !stop_flag)
				  printStatusLine(system_options, current_seed, STOPRESULT_NORMAL, stime, traverse->tag );
				
				stop_flag = true;
			  }
			if( !state_changed && transition_states[idx] != checkresult)
			  state_changed = true;
			transition_states[idx] = checkresult;
			traverse = traverse->next;
		  }
		delete first; // we can do this now as we no longer need to
					  // save the stoplist until the loop exits, due
					  // to moving printStatusLine to immediately upon
					  // finding the stopping condition.
		if( state_changed )
		  {
			sendTransitionStateVectorToPython( transition_states, stime );
			state_changed = false;
		  }
	  }
  } while( stime < maxsimtime && !stop_flag);
  
  if( stime == NAN )
    printStatusLine(system_options, current_seed, STOPRESULT_NAN, 0.0, NULL);
  else if ( stop_flag )
    {
	  dumpCurrentStateToPython();
    }
  else // stime >= maxsimtime
	{
	  dumpCurrentStateToPython();
	  printStatusLine(system_options, current_seed, STOPRESULT_TIME, maxsimtime, NULL );
	}  


  // 	if( stopcount > 0 )
  // 	  {
  // 		curcount = 1;
  // 		checkresult = 0;
  // 		first = getStopComplexList( system_options, 0 );
  // 		checkresult = complexList->checkStopComplexList( first->citem );
  // 		traverse = first;
  // 		while( curcount < stopcount && checkresult == 0 )
  // 		  {
  // 			traverse = traverse->next;
  // 			checkresult = complexList->checkStopComplexList( traverse->citem );
  // 			curcount++;
  // 		  }
  // 		if( checkresult == 0 )
  // 		  {
  // 			delete first;
  // 		  }
  // 	  }
  // 	if( checkresult > 0 )
  // 	  {
  // 		;// printTrajLine(system_options, traverse->tag, stime );
  // 	  }
  // 	else
  // 	  ;
  // 	//          printTrajLine(system_options,"NOSTATE", stime );
  // } while( /*rate > 0.01 && */ stime < maxsimtime && !(checkresult > 0 && stopindex == curcount-1));
  // if( checkresult > 0 )
  // 	delete first;
  // // if( ointerval < 0 || testLongAttr(system_options, output_state ,=, 0 ))
  // //   complexList->printComplexList( 0 );

  // if( stime == NAN )
  // 	printStatusLine(system_options, current_seed, STOPRESULT_NAN, 0.0, NULL);
  // else
  // 	printStatusLine(system_options, current_seed, STOPRESULT_TIME, stime, NULL );
    


}


void SimulationSystem::StartSimulation_FirstStep( void )
{
  InitializeRNG();
  
  while( simulation_count_remaining > 0 )
    {
      if( InitializeSystem() != 0)
		return;

      SimulationLoop_FirstStep();
      generateNextRandom();
      pingAttr( system_options, increment_trajectory_count );
      simulation_count_remaining--;
    }
}



void SimulationSystem::SimulationLoop_FirstStep( void )
{
  double rchoice,rate,stime=0.0, ctime=0.0;
  bool checkresult = false;

  double maxsimtime;
  long stopcount;
  long stopoptions;
  class stopcomplexes *traverse = NULL, *first = NULL;
  long ointerval;
  long trajMode;
  double otime;
  double otime_interval;
  double frate = 0.0;

  getLongAttr(system_options, output_interval,&ointerval);		// FD: here the parameters set, and the PyObject is accessed
  getDoubleAttr(system_options, output_time,&otime);
  getLongAttr(system_options, use_stop_conditions,&stopoptions);
  getLongAttr(system_options, stop_count,&stopcount);
  getDoubleAttr(system_options, simulation_time,&maxsimtime);
  getDoubleAttr(system_options, output_time, &otime_interval );

  complexList->initializeList();

  rate = complexList->getJoinFlux();

  // scomplexlist returns a 0.0 rate if there was a single complex in
  // the system, and a -1.0 rate if there are exactly 0 join moves. So
  // the 0.0 rate should probably be caught, though if you use a
  // single complex system for a starting state it's probably
  // deserved.

  if ( rate == 0.0 )
    { // no initial moves
      printStatusLine_First_Bimolecular( system_options, current_seed, STOPRESULT_NOMOVES, 0.0, 0.0, NULL );
      return;
    }

  rchoice = rate * drand48();

  complexList->doJoinChoice( rchoice );

  // store the forward rate used for the initial step so we can record it.
  frate = rate * dnaEnergyModel->getJoinRate_NoVolumeTerm() / dnaEnergyModel->getJoinRate() ; 

  // Begin normal steps.
  rate = complexList->getTotalFlux();
  do {
    
    rchoice = rate * drand48();
    stime += (log( 1. / (1.0 - drand48()) ) / rate ); 
    

    complexList->doBasicChoice( rchoice, stime );
    rate = complexList->getTotalFlux();

    if( stopcount > 0 && stopoptions)
      {
        checkresult = false;
        first = getStopComplexList( system_options, 0 );
        traverse = first;
        checkresult = complexList->checkStopComplexList( traverse->citem );
        while( traverse->next != NULL && !checkresult )
          {
            traverse = traverse->next;
            checkresult = complexList->checkStopComplexList( traverse->citem );
          }
        if( !checkresult && first != NULL)
          delete first;
      }
  } while( stime < maxsimtime && !checkresult );
      
  // if( !sMode && otime > 0.0 )
  //   {
  //     complexList->printComplexList(0);
  //     printf("Final state reached: Time: %6.6e\n",stime);
  //   }

  if( checkresult )
    {
	  dumpCurrentStateToPython();
      if( strcmp( traverse->tag, "REVERSE") == 0 )
        printStatusLine_First_Bimolecular( system_options, current_seed, STOPRESULT_REVERSE, stime, frate, traverse->tag );
      else 
        printStatusLine_First_Bimolecular( system_options, current_seed, STOPRESULT_FORWARD, stime, frate, traverse->tag );
      delete first;
    }
  else
    {
	  dumpCurrentStateToPython();
      printStatusLine_First_Bimolecular( system_options, current_seed, STOPRESULT_FTIME, stime, frate, NULL );
    }
  // if( ! sMode )
  //   printf("Trajectory Completed\n");
  
}

 
 ///////////////////////////////////////////////////////////
 // void dumpCurrentStateToPython( void );				  //
 // 													  //
 // Helper function to send current state to python side. //
 ///////////////////////////////////////////////////////////

void SimulationSystem::dumpCurrentStateToPython( void )
{
  int id;
  char *names, *sequence, *structure;
  double energy;
  SComplexListEntry *temp;
  temp = complexList->dumpComplexListToPython();
  while( temp != NULL )
	{
	  temp->dumpComplexEntryToPython( &id, &names, &sequence, &structure, &energy );
	  printComplexStateLine( system_options, current_seed, id, names, sequence, structure, energy );
	  temp = temp->next;
	}
}


/////////////////////////////////////////////////////////////////////////////////////
// void sendTransitionStateVectorToPython( boolvector transition_states );		   //
// 																				   //
// Helper function to prepare a Python list object containing the bool information //
//  about which transition states we are in.									   //
/////////////////////////////////////////////////////////////////////////////////////

void SimulationSystem::sendTransitionStateVectorToPython( boolvector transition_states, double current_time)
{
  PyObject *mylist = PyList_New( (Py_ssize_t) transition_states.size() );
  // we now have a new reference here that we'll need to DECREF.

  if( mylist == NULL )
	return; // TODO: Perhaps raise an exception to the Python side here that
			// we couldn't pass the information back...
  
  boolvector_iterator it;
  Py_ssize_t index = 0;

  for( it=transition_states.begin(); it<transition_states.end(); it++ )
	{
	  if( *it ) // bool value was true
		{
		  Py_INCREF( Py_True );
		  PyList_SET_ITEM( mylist, index, Py_True );
		  // ownership of this reference to Py_True has now been stolen by PyList_SET_ITEM.
		}
	  else
		{
		  Py_INCREF( Py_False );
		  PyList_SET_ITEM( mylist, index, Py_False );
		  // ownership of this reference to Py_False has now been stolen by PyList_SET_ITEM.
		}
	  index++;
	}

  PyObject *transition_tuple = Py_BuildValue("dO", current_time, mylist);
  // we now have a new reference to transition_tuple.  note that our
  // reference to mylist has NOT been stolen by this call [it was
  // increffed in the call itself, so we can decref now with no
  // worries]
  Py_DECREF( mylist );
  // we now have no references to mylist directly owned [though transition tuple has one via BuildValue]
  
  pushTransitionInfo( system_options, transition_tuple );
  // transition_tuple has been decreffed by this macro, so we no longer own any references to it

}


 ///////////////////////////////////////////////////////////
 // void sendTrajectory_CurrentStateToPython( void );	  //
 // 													  //
 // Helper function to send current state to python side. //
 ///////////////////////////////////////////////////////////

void SimulationSystem::sendTrajectory_CurrentStateToPython( double current_time )
{
  int id;
  char *names, *sequence, *structure;
  double energy;
  SComplexListEntry *temp;
  temp = complexList->dumpComplexListToPython();
  while( temp != NULL )
	{
	  temp->dumpComplexEntryToPython( &id, &names, &sequence, &structure, &energy );
	  pushTrajectoryComplex( system_options, current_seed, id, names, sequence, structure, energy );
	  temp = temp->next;
	}
  pushTrajectoryInfo( system_options, current_time );
}

///////////////////////////////////////////////////
// This has now been ref count checked, etc etc. //
///////////////////////////////////////////////////

        
int SimulationSystem::InitializeSystem( PyObject *alternate_start )
{
  class StrandComplex *tempcomplex;
  char *sequence, *structure;
  class identlist *id;
  int start_count;
  PyObject *py_start_state = NULL, *py_complex = NULL;
  PyObject *py_seq = NULL, *py_struc = NULL;
  PyObject *py_err = NULL;

  startState = NULL;
  if( complexList != NULL )
    delete complexList;

  complexList = new SComplexList( dnaEnergyModel );

  if( alternate_start != NULL )
    py_start_state = alternate_start;
  else
    py_start_state = getListAttr(system_options, start_state);
  // new reference

  start_count = PyList_GET_SIZE(py_start_state);  
  // doesn't need reference counting for this size call.
  // the getlistattr call we decref later.
  
  for( int index = 0; index < start_count; index++ )
    {
      // #ifndef DEBUG_MACROS
      py_complex = PyList_GET_ITEM(py_start_state, index);
      // Borrowed reference, we do NOT decref it at end of loop.
            
      // #else
      //       py_complex = PyList_GetItem(py_start_state, index);
      // #endif

#ifdef DEBUG_MACROS
      printPyError_withLineNumber();
#endif

      sequence = getStringAttr(py_complex, sequence, py_seq);
      // new reference
      
      structure = getStringAttr(py_complex, structure, py_struc);
      // new reference
	  // Need to check if an error occurred, specifically, it could be an IOError due to sample failing. If so, we need to get the heck out of dodge right now.
	  py_err = PyErr_Occurred(); 
	  // py_err is a borrowed reference
	  if( py_err != NULL )
		{ // then an error occurred while getting the structure. Test for IOError (sample failure):
		  if( PyErr_ExceptionMatches(PyExc_IOError) )
			fprintf(stderr,"MULTISTRAND: Starting Structure could not be retrieved for index %d in your options object's start_state. This is likely due to Boltzmann sampling failing: please check that the program 'sample' exists and points correctly to the NUPACK sample binary. Or try 'print o.start_state[%d].structure' where 'o' is your options object and refer to that error message (if any).\n", index,index );
		  else
			fprintf(stderr,"MULTISTRAND: An unidentified exception occurred while trying to initialize the system.\n");
		  return -1;
		}
	  
      
      id = getID_list( system_options, index, alternate_start );
      
      tempcomplex = new StrandComplex( sequence, structure, id );
      // StrandComplex does make its own copy of the seq/structure, so we can now decref.
      
      Py_DECREF( py_seq );
      Py_DECREF( py_struc );
      startState = tempcomplex;
      complexList->addComplex( tempcomplex );
      tempcomplex = NULL;
    }
  Py_DECREF( py_start_state );
  
  // Update the current seed and store the starting structures
  //   note: only if we actually have a system_options, e.g. no alternate start
  if( alternate_start == NULL && system_options != NULL)
	setLongAttr( system_options, interface_current_seed, current_seed );
  
  return 0;
}

void SimulationSystem::InitializeRNG( void )
{
  bool use_fixed_random_seed = false;
  FILE *fp = NULL;
  getBoolAttr( system_options, initial_seed_flag, &use_fixed_random_seed);

  if( use_fixed_random_seed )
    getLongAttr(system_options, initial_seed,&current_seed);
  else
    {
      if((fp = fopen("/dev/urandom","r")) != NULL )
        {  // if urandom exists, use it to provide a seed
          long deviceseed;
          fread(&deviceseed, sizeof(long), 1, fp);

          current_seed = deviceseed;
          fclose(fp);
        }
      else // use the possibly flawed time as a seed.
        {
          current_seed = time(NULL);
        }
    }
  // now initialize this generator using our random seed, so that we can reproduce as necessary.
  srand48( current_seed );
}

void SimulationSystem::generateNextRandom( void )
{
  current_seed = lrand48();
  srand48( current_seed );
}

PyObject *SimulationSystem::calculateEnergy( PyObject *start_state, int typeflag )
{
  double *values = NULL;
  PyObject *retval = NULL;

  // calc based on current state, do not clean up anything.
  if( start_state != Py_None )
    {
      InitializeSystem( start_state );
      complexList->initializeList();
    }

  values = complexList->getEnergy( typeflag ); // NUPACK energy output : bimolecular penalty, no Volume term.
  // number is complexList->getCount()

  retval = PyTuple_New( complexList->getCount() );
  // New Reference, we return it.
  // The complex list is a linked list and new items are added at the head; so we need to reverse the resulting list to get the data back out.
  for( int loop = complexList->getCount() -1 ; loop >= 0; loop-- )
    PyTuple_SET_ITEM( retval, loop, PyFloat_FromDouble( values[loop] ) );
  // the reference from PyFloat_FromDouble is immediately stolen by PyTuple_SET_ITEM.

  delete[] values;

  return retval;
}
