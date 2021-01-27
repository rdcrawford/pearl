#include "BioSeq.h"
#include <vector>
#include <string>
#include <stdlib.h>

// -----------------------------------------------------------------------------
// Genome
// Ryan D. Crawford
// 05/11/2020
// -----------------------------------------------------------------------------
// This class inherits the propertes of the BioSeq class, representing the
// whole genome sequence. This allows easy parsing of the fasta file to
// retrieve component sequences.
// -----------------------------------------------------------------------------

#ifndef _GENOME_
#define _GENOME_
class Genome: public BioSeq
{
public:

  // Default Ctor
  Genome()
  { ; }

  // Value ctor: takes the paths to the genome annoations file and fasta file
  Genome( const std::string &faPath, const std::string &genomeId ):
    BioSeq( faPath ), genomeId( genomeId )
  { ; }

  // Default Dtor
  ~Genome()
  { ; }

  // Overload the Function to enable sorting of genomes.
  friend bool operator<( const Genome &lhs, const Genome &rhs );

  // Return  the name of this genome
  std::string getGenomeName() const;

  // Return the number of bps in this genome
  unsigned int getGenomeSize() const;


  const Genome* getRef() const;

private:

  // Name of this genome
  std::string genomeId;

};
#endif

// -----------------------------------------------------------------------------
