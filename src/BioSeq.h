#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <iostream>

// -----------------------------------------------------------------------------
// BioSeq
// Ryan D. Crawford
// 05/12/2020
// -----------------------------------------------------------------------------
// This class defines a bioseq class object created by reading in the
// contents of a fasta file as a character vector representing the
// contig, gene sequences, ect. The fasta headers until the first space are
// stored in the "seqNames". If lowercase the sequence is converted to upper
// case.
// -----------------------------------------------------------------------------

#ifndef _BIO_SEQ_
#define _BIO_SEQ_
class BioSeq
{
public:

  // Default ctor
  BioSeq()
  { ; }

  // Value ctor: assignes the path to the fasta file. The file is not read in
  // until requested
  BioSeq( const std::string &faPath ): faPath( faPath )
  { ; }

  // Dtor: does noting
  ~BioSeq()
  { ; }

  // Set the path to the fasta file and parse the data.
  bool setFasta( std::string faPath );

  // This function reads in the fasta file. Returns true if the file
  // was able to be read in
  bool parseFasta();

  // Parses the fasta file from a vector of strings
  bool parseFasta( std::vector< std::string > &faSeq );

  // Write the sequence in a multi-fasta file
  bool writeSeqs( std::string faPath );

  // This function returns the sequence as a vector
  std::vector< std::string > getSeqs();

  // Find the index of the input seq name. If the name is found the "seqIdx"
  // variable is updated and true is returned. False is returned if the
  // seq name is not found
  bool getSeqIndex( const std::string &seqName, unsigned int &seqIdx ) const;

  // Free the memory associated with the whole genome sequence. The contig
  // names are retained. This function exists for instances where
  // it is important to mimize memory usage
  void clearSeqs();

  // Create a vector with the names of the contigs including the fasta headers
  std::vector< std::string > getSeqNames();

  // This function get the substrinconvertToUper();g corresponding to the input
  // parameters the "seq" variable is passed by reference and updated if the
  // parameters are valid. If the parameters are valid and "seq" was sucess
  // updated returns true
  bool getSeqAtCoord( const unsigned int seqIdx, const unsigned int startPos,
    const unsigned int endPos, std::string &seq ) const;

  bool getSeqAtCoord( const std::string &seqName, const unsigned int startPos,
    const unsigned int endPos, std::string &seq ) const;

  // Return the fasta file for these sequences
  std::string getFasta() const;

  // Get the number of contigs in this genome
  int nSeqs() const;

  // Add a sequence to the this BioSeq
  void addSeq( const std::string &faHeader, const std::string &seq );

private:

  // Allow access from genome class and msa class
  friend class Genome;

  // The path to the fasta file corresponding to this genome sequence
  std::string faPath;

  // The sequece contained in the fasta file
  std::vector < std::string > seqs;

  // The names assigned to each contig
  std::vector < std::string > seqNames;

  // Parse the fasta header to get the contig names
  std::string getSeqName( std::string &faHeader );

  // Convert the sequence to upper case
  void convertToUper();

  // Maximium index of the contigs
  unsigned int maxSeqIdx;
};
#endif

// -----------------------------------------------------------------------------
