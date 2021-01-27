#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include <filesystem>
#include <stdlib.h>

// -----------------------------------------------------------------------------
// Input Parser
// Ryan D. Crawford
// 2020/01/15
// -----------------------------------------------------------------------------
// This class parses command line arguemnts and assigns them to the values
// to variables or default values. Then allows access to these variables --
// all are made public.
// -----------------------------------------------------------------------------

class InputParser
{
public:

  // Ctor: takes the command line arguments
  InputParser( int argc, char *argv[] );

  // Vectory with the paths to the fasta files
  std::vector< std::string > fastaFiles;

  // Vector with the genome Ids
  std::vector< std::string > genomeIds;

  // Directory to write temp files
  std::string outDir;

  // Identifer for a run that is appended to the output files
  std::string runId;

  // Output files
  std::string outPath;

  // Minimum identity for
  double minIdent;

  // Minimium lenth of alignment to keep
  unsigned int minLen;

  // Print out the
  void printArgs();

private:

  // Input arguments
  std::vector < std::string > tokens;

  // Look up the option corresponding to a command line argument
  bool getOption( const std::string &option, std::string &param );

  // Update the paths to a directory ( passed by reference ), with the
  // paths to the files in the input directory
  void getPaths( const std::string &dir, std::vector< std::string > &paths,
    const std::string &ext );

  // Find if a parameter exists
  bool findOption( const std::string &option );

  // Print out the command line arguments
  void printOptions();

  // Returns true if a Filesystem::path has the appropriate extension
  inline bool findExt( std::string const &path, std::string const &ext );

  // Extract genome name from path
  std::string getGenomeId( std::string fa, std::string ext );

  // Directory wih the fasta files
  std::string fastaDir;

  // Extension on the fasta files
  std::string fastaExt;
};

// -----------------------------------------------------------------------------
