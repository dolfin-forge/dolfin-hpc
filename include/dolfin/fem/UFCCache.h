// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFC_DATA_H
#define __DOLFIN_UFC_DATA_H

#include <dolfin/common/types.h>
#include <dolfin/fem/UFCCell.h>

namespace dolfin
{

class Form;

/// This class is a simple data structure that holds data used
/// during assembly of a given UFC form. Data is created for each
/// primary argument, that is, v_j for j < r. In addition, nodal
/// basis expansion coefficients and a finite element are created
/// for each coefficient function.

class UFCCache
{
public:
  UFCCache() = default;

  /// Destructor
  ~UFCCache();

  /// initialize datastructure
  void init( Form const & form );

  // Current cell
  UFCCell cell;

  // Local tensor
  std::vector< real > A;

  // Local tensor for macro element
  std::vector< real > macro_A;

  // Array of local dimensions of macro element for primary arguments
  std::vector< size_t > macro_local_dimensions;

  // Array of mapped dofs for primary arguments
  std::vector< size_t * > dofs;

  // Array of mapped dofs of macro element for primary arguments
  std::vector< size_t * > macro_dofs;

  // Array of coefficients
  std::vector< real * > w;

  // Array of coefficients on macro element
  std::vector< real * > macro_w;
};

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_DATA_H */
