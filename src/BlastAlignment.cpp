#include "BlastAlignment.h"

// -----------------------------------------------------------------------------
// BlastAlignment
// Ryan D. Crawford
// 2020/01/13
// -----------------------------------------------------------------------------

// ---- Overload the operators -------------------------------------------------

// Overload the > operator so we can sort blast alignments by position
// in the fasta file
bool operator<( const BlastAlignment &lhs, const BlastAlignment &rhs )
{
  return lhs.length < rhs.length;
}

// Overload the > operator so we can sort blast alignments by the length of
// the alignment
bool operator>( const BlastAlignment &lhs, const BlastAlignment &rhs )
{
  return lhs.length > rhs.length;
}

// Overload the > operator so we can sort blast alignments by position
// in the fasta file
bool operator==( const BlastAlignment &lhs, const BlastAlignment &rhs )
{
  // If these alignments are on different contigs, return which is greater
  // alphabetically
  if ( lhs.qSeqId != rhs.qSeqId ) return false;

  // Return which is further left on the contig
  return lhs.qStart == rhs.qStart && lhs.qEnd == rhs.qEnd;
}

// Overload the != operator so we can sort blast alignments by position
// in the fasta file
bool operator!=( const BlastAlignment &lhs, const BlastAlignment &rhs )
{
  // If these alignments are on different contigs, return which is greater
  // alphabetically
  if ( lhs.qSeqId != rhs.qSeqId ) return true;

  // Return which is further left on the contig
  return lhs.qStart != rhs.qStart || lhs.qEnd != rhs.qEnd;
}

// Overload the + operator so the data on multiple subjects can be stored
// under the same query
void operator+( BlastAlignment &lhs, const BlastAlignment &rhs )
{
  // Append the data to the
  lhs.subjects.insert(
    lhs.subjects.end(), rhs.subjects.begin(), rhs.subjects.end()
    );
}

bool operator==( const Subject & lhs, const Subject & rhs );

// ---- Blast alignment member functions ---------------------------------------

BlastAlignment::BlastAlignment(
  std::stringstream &ss, const Genome* qry, const Genome* subj
  )
{
  // Variable initializtaions for the subject data
  std::string sSeqId;
  int         sStart;
  int         sEnd;
  double      pIdnet;

  // Assign the values from the string strem to the variables
  ss >> qSeqId;
  ss >> sSeqId;
  ss >> qStart;
  ss >> qEnd;
  ss >> sStart;
  ss >> sEnd;
  ss >> length;
  ss >> pIdnet;

  if ( sStart > sEnd )
  {
    auto temp = sEnd;
    sEnd      = sStart;
    sStart    = temp;
  }

  this->query = qry;

  // Create the subject class object
  Subject subject( sSeqId, qStart, qEnd, sStart, sEnd, length, pIdnet, subj );
  subjects.push_back( subject );
}

  // Dtor
BlastAlignment::~BlastAlignment()
{
  // std::cout << "Deleting BlastAlignment!" << std::endl;
  query = nullptr;
}

unsigned int BlastAlignment::getAlgnLen()
{
  return length;
}

void BlastAlignment::printAlign()
{
  for ( int i = 0; i < 80; i++ ) std::cout << '-';
  std::cout << std::endl  << "Query: " << query->getGenomeName() << ": "
            << qSeqId << "; Length: " << length << "; Position: "
            << qStart << "-" << qEnd << std::endl
            << "  -- " << subjects.size() << " subjects "
            << "within this alignment:" << std::endl;
  for ( auto it = subjects.begin(); it != subjects.end(); it++ )
    it->printSubj();
  for ( int i = 0; i < 80; i++ ) std::cout << '-';
  std::cout << std::endl;
}

bool BlastAlignment::checkIsEquiv( const BlastAlignment &qry )
{
  auto it = subjects.begin();
  while ( it != subjects.end() )
  {
    if ( it->sSeqId == qry.qSeqId )
    {
      if ( it->sStart <= qry.qStart && qry.qEnd <= it->sEnd ) return true;
    }
    it++;
  }
  return false;
  // If the subject was not found return false
  // if ( it == qry.subjects.end() ) return false;
  // for ( auto it = qry.subjects.begin(); it != qry.subjects.end(); it++ )
  // {
  //   it = find( subjects.begin(), subjects.end(), *it );
  //   if ( it == subjects.end() ) return false;
  // }
  // return true;
}

bool BlastAlignment::findSubj( const BlastAlignment &qry )
{
  auto it = find( subjects.begin(), subjects.end(), qry.subjects[0] );
  if ( it == subjects.end() ) return false;
  return true;
}

bool  BlastAlignment::getAlignSeq(
  std::string &seqId, std::string &seq
  )
{
  seqId = query->getGenomeName() + ":" + qSeqId + ";" +
    std::to_string( qStart ) + "-" + std::to_string( qEnd );
  return query->getSeqAtCoord( qSeqId, qStart, qEnd, seq );
}

// -----------------------------------------------------------------------------
