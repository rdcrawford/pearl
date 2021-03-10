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
bool operator==( const BlastAlignment &lhs, const BlastAlignment &rhs );


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
  alignments.sort( );
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
    if ( checkIsOverlap( *rIt, *qIt ) )
    {
      // Now we need to handle the subject. This alignment may be
      // within another alignment already present in a subject in this
      // blast alignment class object.
      if ( !rIt->findSubj( *qIt ) )
      {
        *rIt + *qIt;
      } else {

        // Get the next list element
        auto temp = std::next( qIt, 1 );

        // Delete the query alignment
        alignments.erase( qIt );

        // Assign the query to the next alignment
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

bool BlastResults::checkIsOverlap(
  const BlastAlignment &lhs, const BlastAlignment &rhs
  )
{
  // If these alignments are on different contigs, return which is greater
  // alphabetically
  if ( lhs.qSeqId != rhs.qSeqId ) return false;

  // Return if the start position is within the left alignemnt
  if ( lhs.qStart <= rhs.qStart &&  rhs.qStart <= lhs.qEnd ) return true;

  // Return if the end position is within the left alignemnt
  return lhs.qStart <= rhs.qEnd &&  rhs.qEnd <= lhs.qEnd;
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
  for ( int i = 0; i < 40; i++ ) std::cout << "-!";
  std::cout << std::endl;
  algn.printAlign();
  for ( int i = 0; i < 40; i++ ) std::cout << "-!";
  std::cout << std::endl;
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

bool BlastResults::splitAlign(
  std::list< BlastAlignment >::iterator &lhs,
  std::list< BlastAlignment >::iterator &rhs
  )
{
  // Because the left alignment is by definition smaller, and may fall at an
  // equivalent left hand mapping position as the right alignment, we first
  // must acess if the
  if ( rhs->checkWithin( *lhs ) )
  {
    rhs->setStartPos( lhs->qEnd + 1 );
    return true;
  }

  // The alignments are sorted. If the interveining sequnce between
  // the start and end of the alignments is sufficiently long to be it's own
  // sequence, split the alignment at the right subjects end position.
  // return false.
  if ( rhs->qStart - lhs->qStart >= minLen )
  {
    lhs->setEndPos( rhs->qStart );
  }
  return false;
}

void BlastResults::disentangleAlgns()
{
  auto ref    = alignments.begin();
  auto qry    = std::next( ref, 1 );

  // for ( int i = 0; i < 40; i++ ) std::cout << "-|";
  // std::cout << std::endl;
  // Get the unique sequences across all alignments
  int counter = findAlgnSeqs( ref, qry, 0 );

  std::cout << "findAlgnSeqs: " << counter << " iterations complete B-)"
            << std::endl;
  // for ( auto &a : alignments ) a.printAlign();
  // for ( int i = 0; i < 40; i++ ) std::cout << "-|";
  std::cout << std::endl;
}

bool BlastResults::eraseAlgn( std::list< BlastAlignment >::iterator &algn )
{
  auto temp = std::next( algn, 1 );
  alignments.erase( algn );
  algn = temp;
  if ( temp == alignments.end() ) return false;
  return true;
}

bool BlastResults::getNextAlgn( std::list< BlastAlignment >::iterator &algn )
{
  algn = std::next( algn, 1 );
  if ( algn == alignments.end() ) return false;
  return true;
}

bool BlastResults::compareAligns(
  std::list< BlastAlignment >::iterator &ref,
  std::list< BlastAlignment >::iterator &qry
  )
{
  // If these two alignments have no overlap, return false: because alignments
  // are sorted by their right hand mapping position there will be no more
  // alignments to collapse into this one.
  if ( checkIsOverlap( *ref, *qry ) )
  {
    // Check if these two alignments are perfectly equivalent
    if ( *ref == *qry )
    {
      // add this sequence to the firse alignment and get rid of the next
      // alignmentBlastAlignment
      *ref + *qry;
      return true;
    }
    // If the alignment was split and the remaining sequence
    else if ( splitAlign( ref, qry ) )
    {
      // Merge the two alignments
      *ref + *qry;
      return false;
    }
  }
  return false;
}

void BlastResults::getNextRef(
  std::list< BlastAlignment >::iterator &ref,
  std::list< BlastAlignment >::iterator &qry
  )
{
  ref = std::next( ref, 1 );
  while ( ref.length < length ) advance( ref, 1 );
  qry = std::next( ref, 1 );
}

int BlastResults::findAlgnSeqs(
  std::list< BlastAlignment >::iterator ref,
  std::list< BlastAlignment >::iterator qry,
  int counter
  )
{
  // for ( int i = 0; i < 80; i++ ) std::cout << '=';
  // std::cout << std::endl;
  // auto r = std::distance( alignments.begin(), ref );
  // auto q = std::distance( alignments.begin(), qry );
  //
  // std::cout << "findAlgnSeqs function call: " << counter << std::endl
  //           << "  -- ref = " << r << std::endl
  //           << "  -- qry = " << q << std::endl
  //           << "  -- There are "<< alignments.size() << " alignments"
  //           << std::endl;

  counter ++;
  if ( ref == alignments.end() || qry == alignments.end() ) return counter;

  // If the reference and query alignments are equivalent, delete the
  // current query and go on to the next alignment. If the alignmen is split,
  // go to the next alignment. If the query is at the end, update the reference
  // to the next alignment and the query to the subsequent alignment
  if ( compareAligns( ref, qry )  ) eraseAlgn( qry );
  else if ( !getNextAlgn( qry ) ) getNextRef( ref, qry );


  // Before the next function call, make sure the reference alignment and the
  // query have any overlap. If there is no overlap between the two sequences
  // update the reference and query to the next sequence
  if ( !checkIsOverlap( *ref, *qry ) ) getNextRef( ref, qry );

  // Find if there is overlap between the next set
  return findAlgnSeqs( ref, qry, counter ) + counter;
}

// -----------------------------------------------------------------------------
// if ( compareAligns( ref, qry )  )
// {
//   eraseAlgn( qry );
//   q = std::distance( alignments.begin(), qry );
//   std::cout << "  -- Erased: new query " << q << std::endl;
// }
// else if ( !getNextAlgn( qry ) )
// {
//   getNextRef( ref, qry );
//   r = std::distance( alignments.begin(), ref );
//   q = std::distance( alignments.begin(), qry );
//   std::cout << "  -- Updated ref: " << r << std::endl
//             << "  -- Updated qry: " << q << std::endl;
// }
