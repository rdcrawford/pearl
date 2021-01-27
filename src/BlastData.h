#include "BlastResults.h"
#include "BlastAlignment.h"
#include "GenomeData.h"
#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// BlastData
// Ryan D. Crawford
// 2020/01/13
// -----------------------------------------------------------------------------
// This class represents all of the blast results obtained by blasting
// genomes pairwise. The genome data (parsed fasta files) is input to the
// constructor along with a working directory to write temp files and minimum
// identity to collapse alignments. Functionaly is included dereplicate
// the alignments that were generated across multipe genomes and write these
// data to the output file.
// -----------------------------------------------------------------------------

#ifndef _BLAST_DATA_
#define _BLAST_DATA_
class BlastData
{
public:

  // Ctor: take the parsed genome data, and an output directory to write
  // files. Additionally, the minimum identity and minimium length for
  // alignments are input. 
  BlastData( const GenomeData &genomeData, const std::string &outDir,
    const double &minIdent, const unsigned int minLen );

  // Dtor
  BlastData()
  { ; }

  // Get all of the unique sequences in the blast alignments
  void findUniqueAligns( const std::string &outFile );

private:

  // Vector of the blast results for each fasta fle
  std::vector< BlastResults > blastResults;

  // Make a blast database for the input fasta file
  std::string madeBlastDb( const Genome* genome, const std::string &outDir );

  // Use blastn to generate the alignments against the input fasta file and
  // the blast database
  std::string blastFasta( const Genome* query, const Genome* subject,
    const std::string &dbPath, const std::string &outDir );
};
#endif

// -----------------------------------------------------------------------------
