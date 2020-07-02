/***********************************************************************
ParameterReader class is used to simplify the process of reading parameters from an input file and/or from the command line.
Version 1.02 (03-14-2012) Zhi Qiu
***********************************************************************/

#ifndef _ParameterReaderHeader
#define _ParameterReaderHeader

#include <vector>
#include <string>

using std::string;

class ParameterReader
{
  private:
      std::vector<string>* names;
      std::vector<double>* values; // store all parameter names and values
      string removeComments(string str, string commentSymbol); // all substring after "symbol" in "str" will be removed
    void phraseEquationWithoutComments(string equation); // phrase an equation like "x=1", assume string has no comments
    long find(std::string name); // give the index of parameter with "name", or -1 if it does not exist
  public:
    ParameterReader();
    ~ParameterReader();
    void phraseOneLine(string str, string commentSymbol=static_cast<string>("#")); // read and phrase one setting string like "x=1"
    void readFromFile(string filename, string commentSymbol=static_cast<string>("#")); // read in parameters from a file
    void readFromArguments(long argc, char * argv[], string commentSymbol=static_cast<string>("#"), long start_from=1); // read in parameter from argument list. The process starts with index="start_from".
    bool exist(string name); // check if parameter with "name" exists
    void setVal(string name, double value); // set the parameter with "name" to value "value"
    double getVal(string name); // return the value for parameter with "name"
    void echo(); // print out all parameters to the screen
};


#endif

/***********************************************************************
Changelog:
09-20-2011: Ver1.01
 -- Bug fix: If the parameter file that is passed to the readFromFile function does not exist, the program stops instead of going into infinite loops.
03-13-2012: Ver1.02
 -- Bug fix: The commentSymbol parameter was not passed to the phraseOneLine
 function in readFromFile.
***********************************************************************/
