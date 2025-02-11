/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2025 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#ifndef INCLUDE_UTILITY_CC_
#define INCLUDE_UTILITY_CC_

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

#include "basetype.h"
#include "optionlists.h"

using namespace std;


namespace utility {

//helper functions

void sequenceToString(ostream&, BaseType*, int);

char* copyToCharArray(string&);

string copyToString(char*);

void moveType(ostream&, int);

}

#endif /* INCLUDE_UTILITY_CC_ */
