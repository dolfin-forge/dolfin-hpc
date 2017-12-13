// Copyright (C) 2005-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2005-12-19
// Last changed: 2007-04-13

#ifndef __DOLFIN_PARAMETRIZED_H
#define __DOLFIN_PARAMETRIZED_H

#include <dolfin/parameter/Parameter.h>
#include <dolfin/parameter/ParameterList.h>

namespace dolfin
{

  /// This is a base class for parametrized classes. Each object
  /// of a parametrized class holds a local database of parameters
  /// that override the global database of parameters.
  ///
  /// Hierarchies of parametrized objects may be created with the
  /// property that each object inherits the parameters of its parent
  /// and all parameters inherited by the parent. A parameter takes
  /// the value of the first value specified in the hierarchy,
  /// starting at the current object and following the parents back to
  /// the global database.
  ///
  /// Thus, a parameter that has not been specified, using the set()
  /// function, in any given object, takes the value specified for the
  /// parent, whereas any parameter value that has been specified
  /// locally overrides the value specified for the parent.

  class Parametrized
  {

  public:

    /// Constructor
    Parametrized();

    /// Destructor
    virtual ~Parametrized();

    /// Equality
    bool operator==(Parametrized const& other);
    bool operator!=(Parametrized const& other);

    /// Empty
    bool empty() const { return parameters_.empty(); }

    /// Inherit
    bool inherit() const { return (parent_ != NULL); }

    /// Add local parameter
    void add(std::string key, Parameter value);

    /// Set value of local parameter
    void set(std::string key, Parameter value);

    /// Set parent from which to inherit parameters (key must be "parent")
    void set(std::string key, Parametrized const& parent);

    /// Get value of parameter with given key (local or nonlocal)
    Parameter get(std::string key) const;

    /// Check if parameter with given key has been defined locally
    bool has(std::string key) const;

    /// Import parameters from external parameter list
    Parametrized& operator<<(ParameterList const& p);

    /// Export parameters to external parameter list
    Parametrized const& operator>>(ParameterList& p) const;

    /// Import parameters from other
    Parametrized& operator<<(Parametrized const& p);

    /// Export parameters to other
    Parametrized const& operator>>(Parametrized& p) const;

    /// Display parameters
    void disp() const;

  private:

    // Local database of parameters
    ParameterList parameters_;

    // Pointer to parent
    Parametrized const * parent_;

  };

}

#endif
