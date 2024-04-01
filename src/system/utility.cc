/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#include "utility.h"
#include <string>
#include <sstream>
#include "move.h"
#include <energymodel.h>
#include "basetype.h"

char* utility::copyToCharArray(string& myString) {

	char* newArray = (char *) new char[myString.length() + 1];
	strcpy(newArray, myString.c_str());

	return newArray;

}

void utility::sequenceToString(std::ostream& str, BaseType* sequence, int size) {

	// the first and final character is the paired base -- adjust the print for this.

	BaseType preBase = (BaseType) sequence[0];
	BaseType postBase = (BaseType) sequence[size + 1];

	if (preBase < 0 || preBase > 5) {
		cout << "Warning! prebase is outside of range" << endl;
	}
	if (postBase < 0 || postBase > 5) {
		cout << "Warning! postbase is outside of range: " << postBase <<endl;
	}

	str << baseTypeString[(BaseType) sequence[0]] << ":";
	for (int i = 1; i < size + 1; i++)
		str << baseTypeString[(BaseType) sequence[i]];
	str << ":" << baseTypeString[(BaseType) sequence[size + 1]];

}

string utility::copyToString(char* inputCharArray) {

	char* newArray = (char *) new char[strlen(inputCharArray) + 1];
	strcpy(newArray, inputCharArray);

	return string(newArray);

}

void utility::moveType(std::ostream& str, int type) {

	if (type & MOVE_INVALID)
		str << "invalid";
	if (type & MOVE_CREATE)
		str << "create";
	if (type & MOVE_DELETE)
		str << "delete";
	if (type & MOVE_SHIFT)
		str << "shift";
	if (type & MOVE_1)
		str << "_1";
	if (type & MOVE_2)
		str << "_2";
	if (type & MOVE_3)
		str << "_3";

}
