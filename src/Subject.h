#include <string>
#include <iostream>
#include "Genome.h"

// -----------------------------------------------------------------------------
// Subject
// Ryan D. Crawford
// 2020/01/21
// -----------------------------------------------------------------------------
// This cass stores the data on a subject present in the blast alignment.
// -----------------------------------------------------------------------------

#ifndef _SUBJECT_
#define _SUBJECT_
struct Subject
{
  // Default Ctor
  Subject( std::string sSeqId, int qStart, int qEnd, int sStart, int sEnd,
    int length, double pIdent, const Genome* genome ):
    sSeqId( sSeqId ), qStart( qStart ), qEnd( qEnd ), sStart( sStart ),
    sEnd( sEnd ), length( length ), pIdent( pIdent ), genome( genome )
  { ; }

  // Dtor
  ~Subject();

  std::string   sSeqId; // Seq name of the query
  unsigned int  qStart; // Start of alignment in query
  unsigned int  qEnd;   // End of alignment in query
  unsigned int  sStart; // Start of alignment in subject
  unsigned int  sEnd;   // End of alignment in subject
  unsigned int  length; // Alignment length
  double        pIdent; // Percentage of identical matches
  const Genome* genome; // Pointer to the genome of the subject

  void printSubj();
};
#endif

// -----------------------------------------------------------------------------
