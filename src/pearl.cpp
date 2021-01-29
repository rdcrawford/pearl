#include "BlastData.h"
#include "InputParser.h"
#include "GenomeData.h"

// -----------------------------------------------------------------------------
// Pearl
// Ryan D. Crawford
// 2020/01/15
// -----------------------------------------------------------------------------
// This program provides functions to perform blasting of sequences against
// each other. First, the genomic data is parsed and sorted so that that
// genomes of highest quality are handled first. Genomes are then blasted
// against each other. For each query genome, alignments nested within the
// same locus are collapsed for that genome. Finally, the unique alignments
// are identified between all genomes (ie the same sequence is shared by
// more than two genomes). The fasta file with the identified unique sequences
// is then written for down stream analysis.
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
