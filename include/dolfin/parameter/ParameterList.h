// Copyright (C) 2003-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-05-06
// Last changed: 2007-04-13

#ifndef __DOLFIN_PARAMETER_LIST_H
#define __DOLFIN_PARAMETER_LIST_H

#include <dolfin/parameter/Parameter.h>

#include <iomanip>
#include <sstream>

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

  /// Return if empty
  inline bool empty() const { return storage_.empty(); }

  //--- ITERATORS -----------------------------------------------------------
  
  typedef Container::iterator       iterator;
  typedef Container::const_iterator const_iterator;

  inline iterator begin() { return storage_.begin();}
  inline iterator end()   { return storage_.end();}

  inline const_iterator begin() const { return storage_.begin();}
  inline const_iterator end() const   { return storage_.end();}

  // Find key in container
  iterator       find(std::string const& key) { return storage_.find(key); }
  const_iterator find(std::string const& key) const { return storage_.find(key); }

  //-------------------------------------------------------------------------

  /// Equality
  bool operator==(ParameterList const& other) const;
  bool operator!=(ParameterList const& other) const;

  /// Import parameters from other parameter list
  ParameterList& operator<<(ParameterList const& other);

  /// Import parameters from other parameter list
  ParameterList const& operator>>(ParameterList& other) const;

  template<class T>
  static void item(std::string const& k, T v)
  {
    std::stringstream ss;
    ss << std::left << std::setw(32) << k << " = " << v;
    message(ss.str());
  }

  void disp() const
  {
    for (const_iterator it = this->begin(); it != this->end(); ++it)
    {
      ParameterList::item(it->first, it->second);
    }
  }

private:

  // Parameters storage
  Container storage_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_PARAMETER_LIST_H */
