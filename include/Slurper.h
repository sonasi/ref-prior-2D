#ifndef SLURPER_HPP
#define SLURPER_HPP
//-----------------------------------------------------------------------------
//  File:    Slurper.cpp
//  Purpose: Slurp tabular data from a text file
//  Created: 03-May-2005 Harrison B. Prosper
//$Revision: 1.1 $
//-----------------------------------------------------------------------------

#include <fstream>
#include <string>
#include <vector>
#include <map>

#ifdef __WITH_CINT__
#include "TObject.h"
#endif

///
class Slurper
{
 public:

  ///
  Slurper();

  ///
  Slurper(std::string filename, int start=0, int count=0);

  ~Slurper();

  ///
  bool   good();

  ///
  void   close();

  ///
  void   rewind();

  ///
  bool   read();

  /// 
  int    entry();

  ///
  double get(std::string name);

  ///
  void   mget(std::map<std::string, double>& data);

  ///
  int    entries();
  
  ///
  int    size();

  ///
  std::vector<double>      vget(std::vector<std::string>& names);

	///
  std::vector<std::string>      nget();

  ///
  std::vector<std::string> names() {return _name;}

  ///
  bool present(std::string name);

 private:
  std::string _filename;
  int _start;
  int _count;
  int _nrow;
  int _size;
  bool _ok;
  bool _opened;

  std::ifstream*             _stream;
  std::vector<double>        _data;
  std::vector<double>        _buffer;
  std::map<std::string, int> _var;
  std::vector<std::string>   _name;

#ifdef __WITH_CINT__
 public:
  ClassDef(Slurper, 1)
#endif
};

#endif

