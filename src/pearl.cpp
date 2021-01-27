#include "BlastData.h"
#include "InputParser.h"
#include "GenomeData.h"

// -----------------------------------------------------------------------------
//
// Ryan D. Crawford
// 2020/01/15
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

// ---- Main function ----------------------------------------------------------

int main( int argc, char *argv[] )
{

  // Parse the command line arguments and print the inputs
  InputParser inputs( argc, argv );
  inputs.printArgs();

  // Parse the fasta files for the input genomes
  GenomeData genomes( inputs.fastaFiles, inputs.genomeIds );

  // Sort the genomes by the number of contigs and size
  genomes.sortGenomes();

  // Blast the fasta files against each other
  BlastData blastData( genomes, inputs.outDir, inputs.minIdent, inputs.minLen );

  // Blast the genomes and find the alignments that are high identity
  // between sets of genomes
  blastData.findUniqueAligns( inputs.outPath );

  return 0;
}

// -----------------------------------------------------------------------------
