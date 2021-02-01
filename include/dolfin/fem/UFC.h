// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFC_DATA_H
#define __DOLFIN_UFC_DATA_H

#include <dolfin/common/types.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/ufc/ufc.h>

namespace dolfin
{

class Cell;
class DofMapSet;
class Form;
class Mesh;
class MeshDistributedData;

/// This class is a simple data structure that holds data used
/// during assembly of a given UFC form. Data is created for each
/// primary argument, that is, v_j for j < r. In addition, nodal
/// basis expansion coefficients and a finite element are created
/// for each coefficient function.

class UFC
{

public:
  /// Constructor
  UFC( Form const & form );

  /// Destructor
  ~UFC();

  // Form
  ufc::form const & form;

  // Array of finite elements for primary arguments
  std::vector< ufc::finite_element * > finite_elements;

  // Array of finite elements for coefficients
  std::vector< ufc::finite_element * > coefficient_elements;

  // Array of cell integrals
  std::vector< ufc::cell_integral * > cell_integrals;

  // Array of exterior facet integrals
  std::vector< ufc::exterior_facet_integral * > exterior_facet_integrals;

  // Array of interior facet integrals
  std::vector< ufc::interior_facet_integral * > interior_facet_integrals;

  // Current cell
  UFCCell cell;

  // Current pair of cells of macro element
  UFCCell cell0;
  UFCCell cell1;

  // Current pair of local facet indices of macro element
  size_t facet0;
  size_t facet1;

  // Local tensor
  std::vector< real > A;

  // Local tensor for macro element
  std::vector< real > macro_A;

  // Array of local dimensions for each argument
  std::vector< size_t > local_dimensions;

  // Array of local dimensions of macro element for primary arguments
  std::vector< size_t > macro_local_dimensions;

  // Array of local dofmap sizes
  std::vector< size_t > local_sizes;

  // Array of global dimensions for primary arguments
  std::vector< size_t > global_dimensions;

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
