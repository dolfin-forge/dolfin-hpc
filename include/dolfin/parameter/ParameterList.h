// Copyright (C) 2003-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-05-06
// Last changed: 2007-04-13

#ifndef __DOLFIN_PARAMETER_LIST_H
#define __DOLFIN_PARAMETER_LIST_H

#include <dolfin/parameter/Parameter.h>

namespace dolfin
{

/// This class represents a database of parameters, where each
/// parameter is uniquely identified by a string.

class ParameterList
{
  typedef _map<std::string, Parameter> Container;
  typedef std::pair<std::string, Parameter> Item;

public:

  /// Constructor
  ParameterList();

  /// Destructor
  ~ParameterList();

  /// Add parameter
  void add(std::string key, Parameter value);

  /// Set value of parameter
  void set(std::string key, Parameter value);

  /// Get value of parameter with given key
  Parameter get(std::string const& key) const;

  /// Check if parameter with given key has been defined
  bool defined(std::string const& key) const;

  //--- ITERATORS -----------------------------------------------------------
  
  typedef Container::iterator       iterator;
  typedef Container::const_iterator const_iterator;

  inline iterator begin() { return storage_.begin();}
  inline iterator end()   { return storage_.end();}

private:

  // Parameters storage
  Container storage_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_PARAMETER_LIST_H */
