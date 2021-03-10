#include "GenomeData.h"

// -----------------------------------------------------------------------------
// GenomeData
// Ryan D. Crawford
// 05/12/2020
// -----------------------------------------------------------------------------

// ---- Function prototypes ----------------------------------------------------

bool operator<( const Genome &lhs, const Genome &rhs );
bool operator>( const Genome &lhs, const Genome &rhs );

// ---- GenomeData Member functions --------------------------------------------

// Value ctor for inputs of gff files, fasta files, and the correspondingw
// genome names
GenomeData::GenomeData(
  const std::vector< std::string > &faPaths,
  const std::vector< std::string > &genomeIds
  )
{
  // Initialize the vector of genome class objects
  genomeData.reserve( faPaths.size() );
  for ( unsigned int i = 0; i < faPaths.size(); i++ )
    genomeData.push_back( Genome( faPaths[i], genomeIds[i] ) );

  for ( auto it = genomeData.begin(); it != genomeData.end(); it++ )
    it->parseFasta();
}

// Return the gene ids for all of the input genomes
std::vector< std::string > GenomeData::getGenomeIds()
{
  // Initializeth vector to output
  std::vector< std::string > genomeIds;

  // Preallocate suffienct memory
  genomeIds.reserve( genomeData.size() );

  // Get the genome ids and add them to the output vector
  for ( auto it = genomeData.begin(); it != genomeData.end(); it++ )
    genomeIds.push_back( it->getGenomeName() );
  return genomeIds;
}

void GenomeData::sortGenomes()
{
  sort( genomeData.begin(), genomeData.end(), std::greater< Genome >() );

  for ( auto &g : genomeData )
  {
    std::cout << "Genome: " << g.getGenomeName()
              << " " << g.getGenomeSize() << " nts" << std::endl;
  }
}

// Return the gene ids for all of the input genomes
std::vector< std::string > GenomeData::getFaPaths()
{
  // Initializeth vector to output
  std::vector< std::string > fastaPaths;

  // Preallocate suffienct memory
  fastaPaths.reserve( genomeData.size() );

  // Get the paths and add them to the output vector
  for ( auto it = genomeData.begin(); it != genomeData.end(); it++ )
    fastaPaths.push_back( it->getFasta() );

  return fastaPaths;
}

unsigned int GenomeData::getNumGenomes() const
{
  return genomeData.size();
}

const Genome* GenomeData::getGenomeRefAtIdx( const int &idx ) const
{
  return & genomeData[ idx ];
}

// -----------------------------------------------------------------------------
