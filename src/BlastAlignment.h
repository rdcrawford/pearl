#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include "Subject.h"
#include "Genome.h"

// -----------------------------------------------------------------------------
// Blast Alignment
// Ryan D. Crawford
// 2020/01/13
// -----------------------------------------------------------------------------
// This struct stores the data for a blast alignment. Inclued are the positions
// of this alignment in the query, a pointer the the genome from which the
// alignment originates. The subjects contained within this alignment are
// stored in a vector. Subjects may not contain the entire sequce of this
// alignment
// -----------------------------------------------------------------------------

#ifndef _BLAST_ALIGNMENT_
#define _BLAST_ALIGNMENT_
struct BlastAlignment
{
  // Default Ctor
  BlastAlignment()
  { ; }

  // Value Ctor: string stream containing a line in the blast output
  BlastAlignment( std::stringstream &ss, const Genome* qry,
    const Genome* subj );

  // Dtor
  ~BlastAlignment();

  const Genome*          query;    // Pointer to the genome of the query
  std::string            qSeqId;   // Seq name of the subject
  unsigned int           qStart;   // Start of alignment in query
  unsigned int           qEnd;     // End of alignment in query
  unsigned int           length;   // Length of the alignment
  std::vector< Subject > subjects; // Start of alignment in subject

  // Return the length of the alignment
  unsigned int getAlgnLen();

  // Check if the input alignment is contained within this alignment
  bool checkIsEquiv( const BlastAlignment &qry );

  // Print the components of the blast alignment
  void printAlign();

  // Returns true if a sequence from this subject is already represnted
  // in a better alignment in this
  bool findSubj( const BlastAlignment &qry );

  // Return the sequence which was aligned in the query. Create a new seq id
  // for this query sequence. Returns true if the "seq" and "seqId"  variables
  // sucessfully updated.
  bool getAlignSeq( std::string &seqId, std::string &seq );
};
#endif

// -----------------------------------------------------------------------------
