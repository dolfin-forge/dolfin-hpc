// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFC_DATA_H
#define __DOLFIN_UFC_DATA_H

#include <dolfin/common/types.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/fem/UFCMesh.h>

#include <ufc.h>

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
  UFC(Form const& form);

  /// Destructor
  ~UFC();

  // Array of finite elements for primary arguments
  ufc::finite_element** finite_elements;

  // Array of finite elements for coefficients
  ufc::finite_element** coefficient_elements;

  // Array of cell integrals
  ufc::cell_integral** cell_integrals;

  // Array of exterior facet integrals
  ufc::exterior_facet_integral** exterior_facet_integrals;

  // Array of interior facet integrals
  ufc::interior_facet_integral** interior_facet_integrals;

  // Form
  ufc::form const& form;

  // Mesh
  UFCMesh mesh;

  // Current cell
  UFCCell cell;

  // Current pair of cells of macro element
  UFCCell cell0;
  UFCCell cell1;

  // Current pair of local facet indices of macro element
  uint facet0;
  uint facet1;

  // Local tensor
  real* A;

  // Local tensor for macro element
  real* macro_A;

  // Array of local dimensions for each argument
  uint* local_dimensions;

  // Array of local dimensions of macro element for primary arguments
  uint* macro_local_dimensions;

  // Array of local dofmap sizes
  uint* local_sizes;

  // Array of global dimensions for primary arguments
  uint* global_dimensions;

  // Array of mapped dofs for primary arguments
  uint** dofs;

  // Array of mapped dofs of macro element for primary arguments
  uint** macro_dofs;

  // Array of coefficients
  real** w;

  // Array of coefficients on macro element
  real** macro_w;

private:

  void init(ufc::form const& form, Mesh& mesh, DofMapSet const& dof_map_set);

};

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_DATA_H */
