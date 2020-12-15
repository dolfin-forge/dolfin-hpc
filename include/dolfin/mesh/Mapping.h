// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU GPL Version 2.

#ifndef __DOLFIN_MAPPING_H
#define __DOLFIN_MAPPING_H

#include <dolfin/common/types.h>

namespace dolfin
{

class Cell;

/**
 *  DOCUMENTATION:
 *
 *  @class Mapping
 *
 *  @brief Declares an interface for mappings.
 */

class Mapping
{

public:

  /// Update map for current element
  virtual void update(Cell const& cell) = 0;

  /// Map given point from the reference element
  virtual void map_from_reference_cell(real const * xref, real * x) const = 0;

  /// Map given point to the reference element
  virtual void map_to_reference_cell(real const * x, real * xref) const = 0;

protected:

  /// Constructor
  Mapping() = default;

  /// Destructor
  virtual ~Mapping() = default;

};

} /* namespace dolfin */

#endif /* __DOLFIN_MAPPING_H */
