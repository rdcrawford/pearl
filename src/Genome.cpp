#include "Genome.h"

// -----------------------------------------------------------------------------
// Genome
// Ryan D. Crawford
// 05/11/2020
// -----------------------------------------------------------------------------

// ---- Overload operators -----------------------------------------------------

bool operator<( const Genome &lhs, const Genome &rhs )
{
  if ( lhs.nSeqs() != rhs.nSeqs() ) return lhs.nSeqs() < rhs.nSeqs();
  return lhs.getGenomeSize() < rhs.getGenomeSize();
}

bool operator>( const Genome &lhs, const Genome &rhs )
{
  if ( lhs.nSeqs() != rhs.nSeqs() ) return lhs.nSeqs() > rhs.nSeqs();
  return lhs.getGenomeSize() > rhs.getGenomeSize();
}

// ---- Member function definitions --------------------------------------------

unsigned int Genome::getGenomeSize() const
{
  // Iterate over each of the contigs. Get the number of residues and add
  // it to the value to output
  int genomeSize = 0;
  for ( auto it = seqs.begin(); it != seqs.end(); it++ )
    genomeSize += it->size();
  return genomeSize;
}

std::string Genome::getGenomeName() const
{
  return genomeId;
}

const Genome* Genome::getRef() const
{
  return this;
}

// -----------------------------------------------------------------------------
