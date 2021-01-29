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

// ---- Subject member functions -----------------------------------------------

void Subject::printSubj()
{
  std::cout << "     -- "     << genome->getGenomeName() << ": "<< sSeqId
            << "; Position: " << qStart << "-" << qEnd
            << "; Lenght: "   << length << " nts"
            << "; Identity: " << pIdent << "%" << std::endl;
}

Subject::~Subject()
{
  // std::cout << "  -- Deleting subject!" << std::endl;
  genome = nullptr;
}

// -----------------------------------------------------------------------------
