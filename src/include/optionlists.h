/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2025 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#ifndef __OPTIONLISTS_H__
#define __OPTIONLISTS_H__

#define STOPTYPE_STRUCTURE                     0
#define STOPTYPE_BOUND                         1
#define STOPTYPE_DISASSOC                      2
#define STOPTYPE_LOOSE_STRUCTURE               3
#define STOPTYPE_PERCENT_OR_COUNT_STRUCTURE    4

#include <Python.h>

#include <string>

using std::string;

// for PyObject *

// The bodies for these functions are defined in options.cc, but the
// actual classes get used across a bunch of the files so they have been
// separated to not cause huge recompiles.


class strandList {
 public:
  char *id;
  char *seq;
  class strandList *next;
  char *lookup( char *item_id );
  char *name_lookup( char *item_id);
  strandList( char *newid, char *newseq, class strandList *old = NULL );
  ~strandList( void );
};


//FD: Again, a self-rolled linked list.
class identList {
 public:
  // Constructors
  identList( long newuid, const char *newid, class identList *old = NULL);
  ~identList( void );

  // Functions
  std::string toString(void);

  // public variables.
  long uid;
  char *id;
  PyObject *pyo_id;  // needed for correct ref counting dealloc.
  class identList *next;
};

class complexItem {
 public:
  // Constructors
  complexItem(const char *struc, class identList *strands, class complexItem *old = NULL);
  complexItem(const char *struc, class identList *strands, class complexItem *old, int newtype);
  complexItem(const char *struc, class identList *strands, class complexItem *old, int newtype, int newcount);
  ~complexItem(void);

  // functions
  string toString();


  // public variables
  char *structure;
  int type;
  int count; // for use with percentage or count stop types
  class identList *strand_ids;
  class complexItem *next;
};

class stopComplexes {
 public:
  stopComplexes(const char *newtag, class complexItem *newitem, class stopComplexes *old = NULL);
  ~stopComplexes(void);

  char *tag;
  string toString(void);


  class stopComplexes *next;
  class complexItem *citem;
};

#endif
//  __OPTIONLISTS_H__
