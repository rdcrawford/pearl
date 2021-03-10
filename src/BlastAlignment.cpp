#include "BlastAlignment.h"

// -----------------------------------------------------------------------------
// BlastAlignment
// Ryan D. Crawford
// 2020/01/13
// -----------------------------------------------------------------------------

// ---- Function prototypes ----------------------------------------------------

bool operator==( const Subject & lhs, const Subject & rhs );

bool operator<( const Subject & lhs, const Subject & rhs );

bool operator>( const Subject & lhs, const Subject & rhs );

// ---- Overload the operators -------------------------------------------------

// Overload the > operator so we can sort blast alignments by position
// in the fasta file
bool operator<( const BlastAlignment &lhs, const BlastAlignment &rhs )
{
  // If this sequences are on different contigs, return which is less
  // by comparing the two strings
  if ( lhs.qSeqId != rhs.qSeqId ) return lhs.qSeqId < rhs.qSeqId;


  // If these sequences are within each other, return which is smaller
  if ( lhs.qStart == rhs.qStart ) return lhs.length < rhs.length;

  // Return which alignmet is further left on the contig
  return lhs.qStart < rhs.qStart;
}

// Overload the > operator so we can sort blast alignments by the length of
// the alignment
bool operator>( const BlastAlignment &lhs, const BlastAlignment &rhs )
{
  // If this sequences are on different contigs, return which is less
  // by comparing the two strings
  if ( lhs.qSeqId != rhs.qSeqId ) return lhs.qSeqId > rhs.qSeqId;


  // If these sequences are within each other, return which is smaller
  if ( lhs.qStart == rhs.qStart ) return lhs.length > rhs.length;

  // Return which alignmet is further left on the contig
  return lhs.qStart > rhs.qStart;
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

  if ( rhs.qStart < lhs.qStart ) lhs.setStartPos( rhs.qStart );

  if ( rhs.qEnd < lhs.qEnd ) lhs.setStartPos( rhs.qEnd );
}

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

  qStart--;
  qEnd--;
  sStart--;
  sEnd--;

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

bool BlastAlignment::checkIsEquiv( BlastAlignment &qry )
{
  // for ( int i = 0; i < 80; i++ ) std::cout << '=';
  // std::cout << std::endl << " ---- Query: " << std::endl;
  // qry.printAlign();
  // std::cout << " ---- This: " << std::endl;
  // this->printAlign();
  // Initialize an iterator for the subjects in this alignment
  auto it = subjects.begin();

  while ( it != subjects.end() )
  {
    // If this subject is a match for the input query
    if ( it->sSeqId == qry.qSeqId )
    {
      // If the alignmen positions in this subject contain the query alignment
      // positions perfectly, return the true
      if ( it->sStart <= qry.qStart && qry.qEnd <= it->sEnd )
      {
        // std::cout << " ------ Match! " << std::endl;
        return true;
      }
    }
    it++;
  }
  // std::cout << " ------ No Match...." << std::endl;
  // If a perfect alignment was not found, return false
  return false;
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

void BlastAlignment::setStartPos( const unsigned int newStart )
{
  // Set the query end position
  this->qStart = newStart;

  // Calculate the length of the alignment
  calcLen();

  // Reset the start positon in all of the subjects of this alignmente
  // for ( auto & subj : this->subjects ) subj.setStartPos( newStart );
}

void BlastAlignment::setEndPos( const unsigned int newEnd )
{
  // Set the query end position
  this->qEnd = newEnd;

  // Calculate the length of the alignment
  calcLen();

  // Reset the start positon in all of the subjects of this alignment
  // for ( auto & subj : this->subjects ) subj.setEndPos( newEnd );
}

void BlastAlignment::calcLen()
{
  this->length = this->qEnd - this->qStart + 1;
}

bool BlastAlignment::checkWithin( const BlastAlignment &qry )
{
  if ( qry.qStart < this->qStart ) return false;
  if ( qry.qEnd > this->qEnd ) return false;
  return true;
}

// -----------------------------------------------------------------------------
