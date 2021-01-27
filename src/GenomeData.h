#include "Genome.h"
#include <vector>
#include <string>
#include <iostream>

// -----------------------------------------------------------------------------
// GenomeData
// Ryan D. Crawford
// 05/12/2020
// -----------------------------------------------------------------------------
// This class provides a functionality for parsing the data for multiple
// genomes.The data for the parsed genomes is stored in a vector. Additional
// functionality is included for gore the genomes by size and number of contgs
// where the genomes with the fewest contigs and then largest size are
// positioned first in the vector. 
// -----------------------------------------------------------------------------

#ifndef _GENOME_DATA_
#define _GENOME_DATA_
class GenomeData
{
public:

  //
  GenomeData( const std::vector< std::string > &faPaths,
    const std::vector< std::string > &genomeIds );

  // Return the gene ids for all of the input genomes
  std::vector< std::string > getGenomeIds();

  // Return the gene ids for all of the input genomes
  std::vector< std::string > getFaPaths();

  // Sort genomes by number of contigs and genome size
  void sortGenomes();

  // Return the number of genomes contained in this object
  unsigned int getNumGenomes() const;

  // Return the reference of the genome at the desired index in the
  // vector of genomes
  const Genome* getGenomeRefAtIdx( const int &idx ) const;

private:

  // This is a vector of genom class objects
  std::vector< Genome > genomeData;

};
#endif

// -----------------------------------------------------------------------------
