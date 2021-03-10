#include "Subject.h"

// -----------------------------------------------------------------------------
// Subject
// Ryan D. Crawford
// 2020/01/21
// -----------------------------------------------------------------------------

// ---- Overload the operators -------------------------------------------------

bool operator==( const Subject & lhs, const Subject & rhs )
{
  if ( lhs.sSeqId != rhs.sSeqId ) return false;
  if ( lhs.sStart <= rhs.sStart && rhs.sEnd <= lhs.sEnd ) return true;
  return false;
}

bool operator<( const Subject & lhs, const Subject & rhs )
{
  return lhs.qStart < rhs.qStart;
}

bool operator>( const Subject & lhs, const Subject & rhs )
{
  return lhs.qStart > rhs.qStart;
}

// ---- Subject member functions -----------------------------------------------

void Subject::printSubj()
{
  std::cout << "     -- "     << genome->getGenomeName() << ": "
            << sSeqId << " (" << sStart << "-" << sEnd << ")"
            << "; Query Position: " << qStart << "-" << qEnd
            << "; Lenght: "   << length << " nts"
            << "; Identity: " << pIdent << "%" << std::endl;
}

Subject::~Subject()
{
  // std::cout << "  -- Deleting subject!" << std::endl;
  genome = nullptr;
}

void Subject::setStartPos( const unsigned int newStart )
{
  // First we need to know if we are moving the start of the alignment the
  // right or left.
  // If the alignment is moved to the left, subtract the difference from the
  // current start.
  if ( this->qStart > newStart)
  {
    this->sStart = this->sStart - ( this->qStart - newStart );
  }
  // Otherwise this is less than and we need to move the alignment to the
  // right and add the difference to the current start
  else
  {
    this->sStart = this->sStart + ( newStart - this->qStart );
  }
  // Update the current start
  this->qStart = newStart;
  calcLen();
}

void Subject::setEndPos( const unsigned int newEnd )
{
  // First we need to know if we are moving the start of the alignment the
  // right or left.
  // If the alignment is moved to the left, subtract the difference from the
  // current start.
  if ( this->sEnd > newEnd)
  {
    this->sEnd = this->sEnd - ( this->sEnd - newEnd );
  }
  // Otherwise this is less than and we need to move the alignment to the
  // right and add the difference to the current start
  else
  {
    this->sEnd = this->sEnd + ( newEnd - this->sEnd );
  }
  // Update the current start
  this->qEnd = newEnd;
  calcLen();
}

void Subject::calcLen()
{
  this->length =  this->sEnd - this->sStart + 1;
}


// -----------------------------------------------------------------------------
