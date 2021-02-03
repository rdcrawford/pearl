#include "BlastResults.h"

// -----------------------------------------------------------------------------
// BlastData
// Ryan D. Crawford
// 2020/01/13
// -----------------------------------------------------------------------------

// ---- Function prototypes ----------------------------------------------------

// Overload the < and > operators so we can sort blast alignments by the
// length of the alignment
bool operator<( const BlastAlignment &lhs, const BlastAlignment &rhs );
bool operator>( const BlastAlignment &lhs, const BlastAlignment &rhs );
void operator+( BlastAlignment &lhs, const BlastAlignment &rhs );
bool operator==( const Subject & lhs, const Subject & rhs );

// ---- Overload the stream operator -------------------------------------------

bool operator>>( BlastResults &blastResults, BlastAlignment &algn )
{
  // If the end of the alignments has been reached, return false and dont update
  if ( blastResults.it == blastResults.alignments.end() )
  {
    // Reset the position of the iterator, for the next comparison
    blastResults.it = blastResults.alignments.begin();
    return false;
  }
  // Update to the current alignment, and increment the iterator
  algn = *blastResults.it;
  blastResults.it ++;
  return true;
}

// ---- BlastResults member functions ------------------------------------------

BlastResults::BlastResults( const BlastResults &rhs )
{
  // Copy over the vartiables
  alignments = rhs.alignments;
  minIdent   = rhs.minIdent;
  minLen     = rhs.minLen;

  // Set the iterator to the beginning of the alignments
  it = alignments.begin();
}


bool BlastResults::parseBlastData(
  std::string path, const Genome* query, const Genome* subject
  )
{
  std::ifstream ifs;
  std::string   line;

  // Open the input stream and check that it is correct
  ifs.open( path.c_str() );
  if ( !ifs.is_open() || ifs.fail() ) return false;
  if ( ifs.peek() == std::ifstream::traits_type::eof() ) return false;

  // Read in the file line by line
  while ( getline( ifs, line ) )
  {
    // Parse the line in the blast data
    std::stringstream ss( line );
    BlastAlignment    algn( ss, query, subject );

    // If this alignment is sufficiently high identity,
    if ( algn.subjects[0].pIdent >= minIdent ) alignments.push_back( algn );
  }
  return true;
}

void BlastResults::sortBlastResults()
{
  alignments.sort( std::greater< BlastAlignment >() );
  it = alignments.begin();
}

bool BlastResults::collapseNestedAligns()
{
  std::cout << "collapseNestedAligns: " << std::endl
            << "  -- Begining: " << alignments.size() << std::endl;

  // If there are no alignments, or there is only one alignment, return false.
  // this function did not collase any alignments
  if ( alignments.size() <= 1 ) return false;

  // Start the iterators in the
  auto rIt = alignments.begin();
  // auto qIt = rIt + 1;
  auto qIt = std::next( alignments.begin(), 1 );

  // This loop finds if alignmens share the same sequence in the query.
  // First, find if one alignment is completely within another alignment.
  // If they are nested, check the the sequence of the subject. If is NOT
  // If the sequence in this subject is not already represented in this
  // alignment, add it to the current alignment
  do
  {
    int r = std::distance( alignments.begin(), rIt );
    int q = std::distance( alignments.begin(), qIt );
    std::cout << "rIt: " << r << " qIt: " << q << std::endl;

    // Find if the postions of these alignments are completely overlapping
    if ( checkIsNested( *rIt, *qIt ) )
    {
      // Now we need to handle the subject. This alignment may be
      // within another alignment already present in a subject in this
      // blast alignment class object.
      if ( !rIt->findSubj( *qIt ) )
      {
        *rIt + *qIt;
      } else {

        auto temp = std::next( qIt, 1 );


        // Delete the query alignment
        alignments.erase( qIt );

        qIt = temp;
      }
    }
    else // If not nested, go to the next alignment
    {
      if ( qIt != alignments.end() ) qIt ++;
    }

    // If the query is at the end of the vector, update the iterators
    // corresponding to the query and reference
    if ( qIt == alignments.end() )
    {
      if ( rIt != alignments.end() )
      {
        advance( rIt, 1 );
        qIt = std::next( rIt, 1 ); // rIt + 1;
      }
    }
  } while ( rIt != alignments.end() && rIt->length >= minLen );

  if ( rIt != alignments.end() ) alignments.erase( rIt, alignments.end() );

  // alignments.shrink_to_fit();
  this->it = alignments.begin();
  std::cout << "  -- End: " << alignments.size() << std::endl;
  return true;
}

bool BlastResults::checkIsNested(
  const BlastAlignment &lhs, const BlastAlignment &rhs
  )
{
  // If these alignments are on different contigs, return which is greater
  // alphabetically
  if ( lhs.qSeqId != rhs.qSeqId ) return false;

  // Return which is further left on the contig
  return lhs.qStart <= rhs.qStart && rhs.qEnd <= lhs.qEnd;
}

void BlastResults::addAlignment( BlastAlignment algn )
{
  alignments.push_back( algn );
}

bool BlastResults::find( BlastAlignment &algn )
{
  // Set the iterator at the first alignment
  auto algnIt = alignments.begin();

  while ( algnIt != alignments.end() )
  {
    if ( algnIt->checkIsEquiv( algn ) ) return true;
    algnIt ++;
  }
  return false;
}

unsigned int BlastResults::nAligns()
{
  return alignments.size();
}

bool BlastResults::writeAlgnSeqs( const std::string &outFile )
{
  // Initialize a "BioSeq" class object to store the deduplicated sequeces
  BioSeq      alignedSeqs;
  std::string faHeader;
  std::string seq;

  // Iterate over the blast results and extract the aligned sequences
  for ( auto it = alignments.begin(); it != alignments.end(); it++ )
  {
    it->getAlignSeq( faHeader, seq );
    alignedSeqs.addSeq( faHeader, seq );
  }

  // Write the fasta file with the unique sequences
  return alignedSeqs.writeSeqs( outFile );
}

// -----------------------------------------------------------------------------
// {
//   for ( int i = 0; i < 80; i++ ) std::cout << '=';
//   std::cout << std::endl;
//   std::cout << "*** Reference alignment ***" << std::endl;
//   algnIt->printAlign();
//   std::cout << "*** Equivalent alignment ***" << std::endl;
//   algn.printAlign();
//   for ( int i = 0; i < 80; i++ ) std::cout << '=';
//   std::cout << std::endl;
//
//   return true;
// }
// bool BlastResults::find( BlastAlignment &algn )
// {
//   // Set the iterator at the first alignment
//   auto algnIt = alignments.begin();
//
//   while ( algnIt != alignments.end() )
//   {
//     if ( algn.length > algnIt->length )
//     {
//       if ( algn.checkIsEquiv( *algnIt ) )
//       {
//         algn.checkIsEquiv( *algnIt );
//         *algnIt = algn;
//         return true;
//       }
//
//     } else {
//
//       if ( algnIt->checkIsEquiv( algn ) ) return true;
//     }
//
//     algnIt ++;
//   }
//   return false;
// }
