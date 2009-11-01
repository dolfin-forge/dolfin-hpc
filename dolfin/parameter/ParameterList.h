// Copyright (C) 2003-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-05-06
// Last changed: 2007-04-13

#ifndef __PARAMETER_LIST_H
#define __PARAMETER_LIST_H

#include <dolfin/common/types.h>
#include "Parameter.h"

namespace dolfin
{
  
  /// This class represents a database of parameters, where each
  /// parameter is uniquely identified by a string.
  
  class ParameterList
  {
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
    Parameter get(std::string key) const;

    /// Check if parameter with given key has been defined
    bool defined(std::string key) const;

    /// Friends
    friend class XMLFile;
    
  private:

    // Parameters stored as an STL map
    _map<std::string, Parameter> parameters;

    // Typedef of iterators for convenience
    typedef _map<std::string, Parameter>::iterator iterator;
    typedef _map<std::string, Parameter>::const_iterator const_iterator;
    
    // Typedef of pair for convenience
    typedef std::pair<std::string, Parameter> pair;
    
  };
  
}

#endif
