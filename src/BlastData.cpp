#include "BlastData.h"

// -----------------------------------------------------------------------------
// Blast Data
// Ryan D. Crawford
// 2020/01/13
// -----------------------------------------------------------------------------

// ---- Function prototypes ----------------------------------------------------

// Overloaded operators used for alignment comparisons:
bool operator<( const BlastAlignment &lhs, const BlastAlignment &rhs );
bool operator>( const BlastAlignment &lhs, const BlastAlignment &rhs );

// Overload the stream operator to iterate over alignments
bool operator >>( BlastResults &blastResults, BlastAlignment &algn );

// ---- BlastData member functions ---------------------------------------------

BlastData::BlastData(
  const GenomeData &genomeData,
  const std::string &outDir,
  const double &minIdent,
  const unsigned int minLen
  )
{
  // Get the number of genomes in this dataset
  unsigned int nGenomes = genomeData.getNumGenomes();

  // Make blast databases
  std::vector< std::string > blastDbs( nGenomes );
  for ( unsigned int i = 1; i < nGenomes; i++ )
    blastDbs[i] = madeBlastDb( genomeData.getGenomeRefAtIdx( i ), outDir );

  // Blast each pair of genomes.
  blastResults.reserve( nGenomes - 1 );
  const Genome* query;
  const Genome* subject;
  for ( unsigned int i = 0; i < nGenomes - 1; i ++ )
  {
    // Create a "blast results" class object and add it to the vector
    blastResults.push_back( BlastResults( minIdent, minLen ) );

    // Get the reference to the query
    query = genomeData.getGenomeRefAtIdx( i );

    // Blast the next set of genomes
    for ( unsigned int j = i + 1; j < nGenomes; j++ )
    {
      subject =  genomeData.getGenomeRefAtIdx( j );
      std::string tsv =
        blastFasta( query, subject, blastDbs[j], outDir );
      blastResults[i].parseBlastData( tsv, query, subject );
    }
  }

  // Sort the blast results by the position in the query
  for ( auto &results : blastResults ) results.sortBlastResults();

  // Find all alignments that are perfectly within another alignment
  for ( auto &results : blastResults ) results.collapseNestedAligns();

  // Set the pointers to null before they go out of scope
  query   = nullptr;
  subject = nullptr;
}


// Create the blast command to run
std::string BlastData::blastFasta(
  const Genome* query,
  const Genome* subject,
  const std::string &dbPath,
  const std::string &outDir
  )
{
  std::string outTsv =
    outDir + query->getGenomeName() + "_" + subject->getGenomeName() + ".tsv";
  std::string blastCmd = "blastn -query " + query->getFasta() +
    " -db " + dbPath + " -out " + outTsv +
     " -outfmt \"6 qseqid sseqid qstart qend sstart send length pident\"";
  system( blastCmd.c_str() );
  return outTsv;
}

// Make Blast DB
std::string BlastData::madeBlastDb(
  const Genome* genome, const std::string &outDir
  )
{
  // Create the path for the blast database to create
  std::string dbPath = outDir + genome->getGenomeName();

  // Make the blast database
  std::string cmd =
    "makeblastdb -dbtype nucl -in " + genome->getFasta() + " -out " + dbPath;
  system( cmd.c_str() );

  return dbPath;
}

void BlastData::findUniqueAligns( const std::string &outFile )
{
  // Initialize to the unique alignments to the first set of alignments is
  // these data. Copy constructore makes a deep copy.
  BlastResults   uniqueAlgns = blastResults[0];
  BlastAlignment algn;

  std::cout << "findUniqueAligns: " << std::endl
            << "  -- Begin: " << uniqueAlgns.nAligns() << std::endl;
  for ( unsigned int i = 1; i < blastResults.size(); i++ )
  {
    while ( blastResults[i] >> algn )
    {
      if ( !uniqueAlgns.find( algn ) ) uniqueAlgns.addAlignment( algn );
    }
  }
  uniqueAlgns.sortBlastResults();
  std::cout << "findUniqueAligns: " << std::endl
            << "  -- End: " << uniqueAlgns.nAligns() << std::endl;

  // Write a fasta file for these sequences
  uniqueAlgns.writeAlgnSeqs( outFile + "_seqs.fasta" );
}

// -----------------------------------------------------------------------------
