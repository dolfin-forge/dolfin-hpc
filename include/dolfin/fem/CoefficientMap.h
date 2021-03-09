// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_COEFFICIENT_MAP_H
#define __DOLFIN_COEFFICIENT_MAP_H

#include <dolfin/fem/Coefficient.h>

#include <string>

namespace dolfin
{

class CoefficientMap : public _ordered_map<std::string, Coefficient *>
{
public:
  /// Check if coefficient label is in the map
  auto has(std::string const& label) const -> bool;

  /// Display information
  void disp() const;

};

}

#endif /* __DOLFIN_COEFFICIENT_MAP_H */
