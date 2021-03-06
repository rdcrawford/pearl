#include "BioSeq.h"
using namespace std;

// -----------------------------------------------------------------------------
// GenomeFeatures
// Ryan D. Crawford
// 05/12/2020
// -----------------------------------------------------------------------------

// This function reads in the fasta file. Returns true if the file
// was able to be read in
bool BioSeq::parseFasta( )
{
  std::string line;        // Current line in the fasta file
  std::string *curContig;  // Pointer to the current contig
  ifstream    ifs;         // File stream for the fasta file
  int         seqIdx = -1; // Index of the contig being parsed

  // If the file stream has failed to open, return false
  ifs.open( faPath.c_str() );
  if ( ifs.fail() || !ifs.is_open() )
  {
     std::cout << "Failed to open the fasta file..." << std::endl;
     exit(1);
  }

  if ( ifs.peek() == std::ifstream::traits_type::eof() ) return false;

  // Read in the fasta file line by line
  while ( getline( ifs, line ) )
  {

    if ( line[0] == '>')
    {
      seqIdx ++; // Advance the couter to the current contig
      // Parse the line to get the name of the contig
      seqNames.push_back( getSeqName( line ) );

      // Create an empty contig to add sequences to
      seqs.push_back( "" );
      curContig = & seqs[ seqIdx ];

    } else {

      // Add the line to the current contig
      *curContig += line;
    }
  }
  // Close up the output file stream
  ifs.close();

  // Set the number of sequences included in the fasta file
  maxSeqIdx = seqIdx;

  // Convert the sequence to upper-case
  convertToUper();

  // Finished. Return true to indicate the the genome was parsed properly
  return true;
}

// This function reads in the fasta file. Returns true if the file
// was able to be read in
bool BioSeq::parseFasta( std::vector< std::string > &faSeq )
{
  std::string *curContig;           // Pointer to the current contig
  ifstream    ifs;                  // File stream for the fasta file
  auto        line = faSeq.begin(); // Iterator for the line in the fasta file
  int         seqIdx = -1;          // Index of the contig being parsed

  // If the file stream has failed to open, return false
  ifs.open( faPath.c_str() );
  if ( ifs.fail() ) return false;

  // Read in the fasta file line by line
  while ( line < faSeq.end() )
  {
    // Check if the first position of the string is a new header
    if ( line->at( 0 ) == '>')
    {
      seqIdx ++; // Advance the couter to the current contig

      // Parse the line to get the name of the contig
      seqNames.push_back( getSeqName( *line ) );

      // Create an empty contig to add sequences to
      seqs.push_back( "" );
      curContig = & seqs[ seqIdx ];
    } else {

      // Add the line to the current contig
      *curContig += *line;
    }
  }
  ifs.close();

  // Set the number of sequences included in the fasta file
  maxSeqIdx = seqIdx;

  // Convert the sequence to upper-case
  convertToUper();

  // Finished. Return true to indicate the the genome was parsed properly
  return true;
}

// Create a vector with the names of the contigs including the fasta headers
std::vector< std::string > BioSeq::getSeqNames()
{
  return seqNames;
}

void BioSeq::convertToUper()
{
  // Any down stream steps with this structure will require it to be in
  // upper case. Test if the first character is lowercase and if it is,
  // transform the entire string to uppercase
  if ( std::islower( seqs[0][0] ) )
  {
    for ( auto &seq : seqs )
      transform( seq.begin(), seq.end(), seq.begin(), ::toupper );
  }
}

// Set the path to the fasta file and parse the data.
bool BioSeq::setFasta( std::string faPath )
{
  this->faPath = faPath;
  return parseFasta();
}

std::string BioSeq::getSeqName( std::string &faHeader )
{
  // Erase the fasta header.
  faHeader.erase( 0, 1 );

  // If the fasta header starts with accn, delete it
  auto accnPos = faHeader.find( "accn|" );
  if ( accnPos != string::npos ) faHeader.erase( 0, accnPos + 5 );

  // Find the position of the first space.
  auto spacePos = faHeader.find( ' ' );

  // If there are no spaces, return the header as is
  if ( spacePos == string::npos ) return faHeader;

  // Erase everything after the space from the header and return
  return faHeader.substr( 0, spacePos);
}

// This function returns the whole genome sequence as a vector
std::vector< std::string > BioSeq::getSeqs( )
{
  return seqs;
}


// Find which contig
bool BioSeq::getSeqIndex(
  const std::string &seqName, unsigned int &seqIdx
  ) const
{
  auto it = find( seqNames.begin(), seqNames.end(), seqName );
  // If the name was not found return false
  if ( it == seqNames.end() ) return false;
  seqIdx = std::distance( seqNames.begin(), it );
  return true;
}

// Free the memory associated with the whole genome sequence. The contig
// names are retained. This function exists for instances where mimizing memory
// usage is critical and keeping the whole genome sequence in memory
void BioSeq::clearSeqs()
{
  seqs.clear();
}

// This function parses sequences to retrieve the substring corresponding
// to the input coordinates. The "seq" variable is passed by reference and
// updated if the input coordinates are valid. If "seq" was successfully
// updated returns true, false otherwise.
bool BioSeq::getSeqAtCoord(
  const unsigned int seqIdx, const unsigned int startPos,
  const unsigned int endPos, std::string &seq
  ) const
{
  // Check that this index corresponds to a contig in this file
  if ( seqIdx > maxSeqIdx ) return false;

  // Check that the requested direction is correct
  if ( startPos >= endPos ) return false;

  // Check the the end of the sequence is not out of range
  if ( endPos > seqs[ seqIdx ].size() ) return false;

  // Update the sequence to return the substring for this contig. The start
  // is decremented by one because it is one indexed
  int len = endPos - startPos + 1;

  // Subset the
  seq = seqs[ seqIdx ].substr( startPos, len );

  return true;
}

bool BioSeq::getSeqAtCoord(
  const std::string &seqName, const unsigned int startPos,
    const unsigned int endPos, std::string &seq
  ) const
{
  // Look up the name of the contig
  unsigned int seqIdx;
  if ( !getSeqIndex( seqName, seqIdx ) ) return false;

  // Check that the requested direction is correct
  if ( startPos >= endPos ) return false;

  // Check the the end of the sequence is not out of range
  if ( endPos > seqs[ seqIdx ].size() ) return false;

  // Update the sequence to return the substring for this contig. The start
  // is decremented by one because it is one indexed
  int len = endPos - startPos + 1;

  // Subset the
  seq = seqs[ seqIdx ].substr( startPos, len );

  return true;
}

void BioSeq::addSeq( const std::string &faHeader, const std::string &seq )
{
  seqs.push_back( seq );
  seqNames.push_back( faHeader );
}

// Write the sequence in a multi-fasta file
bool BioSeq::writeSeqs( std::string faPath )
{
  // Initialize the output file stream and open for writing
  std::ofstream ofs;
  ofs.open( faPath.c_str() );
  if ( ofs.fail() || !ofs.is_open() ) return false;

  for ( unsigned int i = 0; i < seqs.size(); i++ )
    ofs << ">" << seqNames[ i ] << endl << seqs[ i ] << endl;

  ofs.close();
  return true;
}

std::string BioSeq::getFasta() const
{
  return faPath;
}

int BioSeq::nSeqs() const
{
  return maxSeqIdx;
}

// -----------------------------------------------------------------------------
