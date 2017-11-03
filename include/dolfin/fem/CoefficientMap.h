// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-11-10
// Last changed: 2014-11-10

#ifndef __DOLFIN_COEFFICIENT_MAP_H
#define __DOLFIN_COEFFICIENT_MAP_H

#include <dolfin/common/Array.h>
#include <dolfin/fem/Coefficient.h>

#include <map>
#include <string>

namespace dolfin
{

class CoefficientMap
{
  typedef std::map<std::string, Coefficient *>  Container;
  typedef std::pair<std::string, Coefficient *> Item;

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

  /// Accessor
  inline Coefficient * operator[](std::string const& label)
  {
    return map_[label];
  }

  /// Accessor (const)
  inline Coefficient * operator[](std::string const& label) const
  {
    Container::const_iterator it = map_.find(label);
    if (it == map_.end()) return NULL;
    return it->second;
  }

  /// Return the size of the coefficient map
  uint size() const;

  /// Clear coefficient map
  void clear();

  /// Display information
  void disp() const;

private:

  Container map_;

};

}

#endif /* __DOLFIN_COEFFICIENT_MAP_H */
