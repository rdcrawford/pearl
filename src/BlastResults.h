#include "BlastAlignment.h"
#include "Genome.h"
#include "BioSeq.h"
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <string>

// -----------------------------------------------------------------------------
// blastResults
// Ryan D. Crawford
// 2020/01/13
// -----------------------------------------------------------------------------
// This class stores "BlastAlignment" class objects in  a vector. The
// blast results in custom output format 6 are parsed and stored in the
// vector.
// -----------------------------------------------------------------------------

#ifndef _BLAST_RESULTS_
#define _BLAST_RESULTS_
class BlastResults
{
public:

  // Ctor: takes the path to the blast resutls in outformat 6 and a
  BlastResults()
  { ; }

  // Ctor: takes the path to the blast resutls in outformat 6 and a
  BlastResults( double minIdent, unsigned int minLen ):
    minIdent( minIdent ), minLen( minLen )
  { ; }

  // Dtor
  ~BlastResults()
  { ; }

  // Copy constructor
  BlastResults( const BlastResults &rhs );

  // Sort the vector of alignments by the position in the fasta file
  void sortBlastResults();

  // Read in the tsv and parse the data into a vector of alignment Class
  // objects
  bool parseBlastData( std::string path, const Genome* query,
    const Genome* subject );

  // Identify alignments that fall within each other and
  bool collapseNestedAligns();

  // Add an alignment to the vector of alignments
  void addAlignment( BlastAlignment algn );

  // Search over the blast alignments and determine if there is an equivalent
  // alignment present in these alignments
  bool find( BlastAlignment &algn );

  // Return the number of alignments in
  unsigned int nAligns();

  // Stream operator for accessing alignments
  friend bool operator>>( BlastResults &blastResults, BlastAlignment &algn );

  // Write a multi fasta file with the aligned sequeces
  bool writeAlgnSeqs( const std::string &outFile );

private:

  // Vector with each individual
  std::list< BlastAlignment > alignments;

  // Iterator to keep track of the current alignment
  std::list< BlastAlignment >::iterator it = alignments.begin();

  // Minimium identity threshold to keep the alignment
  double minIdent;

  // Minimium lenth for the query alignment
  unsigned int minLen;

  // Function to identify an alignment that is with another alignment
  bool checkIsNested( const BlastAlignment &lhs, const BlastAlignment &rhs );

};
#endif

// -----------------------------------------------------------------------------
