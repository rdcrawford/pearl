#include "InputParser.h"
namespace fs = std::filesystem;
using namespace std;

// -----------------------------------------------------------------------------
//
// Ryan D. Crawford
// 2020/01/15
// -----------------------------------------------------------------------------

InputParser::InputParser ( int argc, char *argv[] )
{
  // Generate a string vector for the command line arguments
  for ( int i = 1; i < argc; i++ )
    this->tokens.push_back( string( argv[i] ) );

  if ( argc < 2 || findOption( "--help" ) )
  {
    printOptions();
    exit( 0 );
  }

  // Get the directory with the fasta files
  if ( !getOption( "--fastaDir", fastaDir ) )
  {
    cout << "Missing reqired argument: --fastaDir" << endl;
    exit( 1 );
  }

  // Get the directory with the fasta files
  if ( !getOption( "--fastaExt", fastaExt ) ) fastaExt = ".fasta";

  // Get the paths to the input files
  getPaths( fastaDir, fastaFiles, fastaExt );

  // Set the genome Ids
  genomeIds.reserve( fastaFiles.size() );
  for ( const auto & fa : fastaFiles )
    genomeIds.push_back( getGenomeId( fa, fastaExt ) );

  if ( fastaFiles.size() <=1 )
  {
    cout << fastaFiles.size()
         << " genomes were input. At least two files must be input"
         << endl;
    exit( 1 );
  }

  // Parse the optional arguments
  string val;
  if ( getOption( "--minIdent", val ) )
  {
    minIdent = std::stod( val );

    if ( minIdent > 1 || minIdent < 0 )
    {
      cout << "Argument --minIdent must be in the range: 0 and 1"
           << endl;
      exit( 1 );
    }
  }
  else // Assign to the minimium
  {
    minIdent = 99;
  }

  if ( !getOption( "--outDir", outDir ) )
  {
    outDir = "";
  } else {

    if ( outDir.back() != '/' ) outDir = outDir + '/';
  }

  if ( !getOption( "--runId", runId ) )
  {
    runId = "";
  }
  outPath = outDir + "pearl_" + runId;


   if ( !getOption( "--minLen", val ) )
   {
     minLen = 500;
   } else {
     minLen = stoi( val );
   }
}

void InputParser::printOptions()
{
  cout << endl << "Required arguments:" << endl <<
     // fastaDir
     "  --fastaDir  Directory containing fasta files to compare"  << endl <<
     endl << "Optional arguments:" << endl
     << "  --minIdent Minimium identity to collapse alignment to "
     << "(range 0-100) Defaults to 99%" << endl
     << "  --outDir   Directory to write files."
     << "defaults to current working directory"  << endl
     << "  --fastaExt Extension for the fasta files. Defaults to .fasta" << endl
     << "  --runId    String to prepend to the output files"  << endl
     << "  --minLen   Minimium length of the alignments to keep. Defaults to "
     << "500 nts."
     << endl << endl;
}

void InputParser::printArgs()
{
  cout << endl << "Pearl:" << endl
    << "  -- " << "Finding conserved alignments for at an identity of "
    << minIdent << "% and at least " << minLen << " nts "<< fastaFiles.size()
    << " genomes" << endl;

  if ( outDir == "" )
  {
    cout << "  -- Writing results to the current working directory" << endl;
  }
  else
  {
    cout << "  -- Writing results to " << outDir << endl;
  }
}

bool InputParser::getOption( const string &option, string &param )
{
  auto it = find( this->tokens.begin(), this->tokens.end(), option );
  if ( it == this->tokens.end() ) return false;
  it++;
  param = *it;
  return true;
}

bool InputParser::findOption( const string &option )
{
  auto it = find( this->tokens.begin(), this->tokens.end(), option );
  if ( it == this->tokens.end() ) return false;
  return true;
}

inline bool InputParser::findExt(
  std::string const &path, std::string const &ext
  )
{
  if ( ext.size() > path.size() ) return false;
  return std::equal( ext.rbegin(), ext.rend(), path.rbegin() );
}

void InputParser::getPaths(
  const string &dir, vector< string > &paths, const std::string &ext
  )
{
  fs::recursive_directory_iterator it = fs::recursive_directory_iterator( dir );
  for ( const auto & dirEntry : it )
  {
    if ( findExt( dirEntry.path(), ext ) ) paths.push_back( dirEntry.path() );
  }
}

string InputParser::getGenomeId( string fa, string ext )
{
  int start = fa.find_last_of( "/" ) + 1;
  int len   = fa.find( ext ) - start;
  return fa.substr( start, len );
}

// -----------------------------------------------------------------------------
