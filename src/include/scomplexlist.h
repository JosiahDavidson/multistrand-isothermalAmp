/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2024 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#ifndef __SCOMPLEXLIST_H__
#define __SCOMPLEXLIST_H__

#include <stdio.h>
#include <iostream>

#include "scomplex.h"

using std::cout;

class SComplexListEntry;
class JoinCriterea;

class SComplexList {
public:

	SComplexList(EnergyModel *energyModel);
	~SComplexList(void);

	SComplexListEntry *addComplex(StrandComplex *newComplex);
	void initializeList(void);
	void regenerateMoves(EnergyModel *energyModel);
	int getMoveCount(void);
	double getMoveFlux(void);
	double getJoinFlux(void);
	double getTotalFlux(void);

	BaseCount getExposedBases();
	OpenInfo getOpenInfo();

	double getJoinFluxArr(void);
	double computeArrBiRate(SComplexListEntry*);
	double cycleCrossRateArr(StrandOrdering*, StrandOrdering*);
	int getCount(void);
	double *getEnergy(int volume_flag);

	void printComplexList(std::ostream&);
	void printStructures(std::ostream&);
	void printAllMoves(std::vector<string>&);
	SComplexListEntry *getFirst(void);
	double doBasicChoice(
		SimTimer&, int* = NULL, int* = NULL, int* = NULL,
		std::stringstream* = NULL, std::stringstream* = NULL);
	JoinCriteria cycleForJoinChoice(SimTimer&, int* = NULL, int* = NULL);
	JoinCriteria cycleForJoinChoiceArr(SimTimer&, int* = NULL, int* = NULL);
	JoinCriteria findJoinNucleotides(
		BaseType, int, BaseCount&, SComplexListEntry*, int* = NULL, HalfContext* = NULL);
	double doJoinChoice(
		SimTimer&, int* = NULL, int* = NULL,
		std::stringstream* = NULL, std::stringstream* = NULL);
	bool checkStopComplexList(class complexItem *stoplist);
	void toString(std::ostream&);
	void updateOpenInfo(void);

private:
	bool checkStopComplexList_Bound(class complexItem *stoplist);
	bool checkStopComplexList_Structure_Disassoc(class complexItem *stoplist);
	bool checkLooseStructure(const char *our_struc, const char *stop_struc, int count);
	bool checkCountStructure(const char *our_struc, const char *stop_struc, int count);

	int numOfComplexes = 0;
	int idcounter = 0;

	SComplexListEntry* first = NULL;
	EnergyModel* eModel = NULL;

	double joinRate = 0.0;	// joinrate is the sum of collision rates in the state

};

class SComplexListEntry {
public:
	SComplexListEntry(StrandComplex *newComplex, int newid);
	~SComplexListEntry(void);
	void initializeComplex(EnergyModel *energyModel);
	void regenerateMoves(EnergyModel *energyModel);
	void fillData(EnergyModel *em);
	void toString(std::ostream&, EnergyModel *);
	void dumpComplexEntryToPython(ExportData& data, EnergyModel *energyModel);

	int id;
	StrandComplex* thisComplex;
	energyS ee_energy;
	double energy;
	double rate;

	SComplexListEntry *next;
};

#endif
