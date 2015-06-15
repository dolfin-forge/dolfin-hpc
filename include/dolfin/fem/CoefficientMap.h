// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-11-10
// Last changed: 2014-11-10

#ifndef __COEFFICIENT_MAP_H_
#define __COEFFICIENT_MAP_H_

#include <dolfin/common/Array.h>
#include <dolfin/function/Function.h>

#include <map>
#include <string>

namespace dolfin
{

class CoefficientMap
{

public:

  /// Constructor
  CoefficientMap();

  /// Destructor
  ~CoefficientMap();

  /// Constructor
  CoefficientMap(CoefficientMap const& other);

  /// Assignment
  CoefficientMap& operator=(CoefficientMap const& other);

  /// Check if coefficient label is in the map
  bool has(std::string const& label) const;

  /// Get coefficient function mapped from the given label
  Function * get(std::string const& label) const;

  /// Return the size of the coefficient map
  uint size() const;

  /// Set coefficient label to map to the given function
  void set(std::string const& label, dolfin::Function& coefficient);

  /// Clear coefficient map
  void clear();

  /// Display information
  void disp() const;

private:

  std::map<std::string, dolfin::Function *> map_;

};

}

#endif /* __COEFFICIENT_MAP_H_ */
