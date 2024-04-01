/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#ifndef INCLUDE_UTILITY_CC_
#define INCLUDE_UTILITY_CC_

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

#include "optionlists.h"
#include "basetype.h"

using namespace std;


namespace utility {

//helper functions

void sequenceToString(std::ostream&, BaseType*, int);

char* copyToCharArray(string&);

string copyToString(char*);

void moveType(std::ostream&, int);

}

#endif /* INCLUDE_UTILITY_CC_ */
