// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-01-17
// Last changed: 2008-06-10

#ifndef __UFC_DATA_H
#define __UFC_DATA_H

#include "UFCMesh.h"
#include "UFCCell.h"

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

  /// Constructor
  UFC(ufc::form const& form, Mesh& mesh, DofMapSet const& dof_map_set);

  /// Destructor
  ~UFC();

  /// Update current cell
  void update(Cell& cell, MeshDistributedData& distdata);

  /// Update current pair of cells for macro element, global numbering
  void update(Cell& cell0, Cell& cell1, MeshDistributedData& distdata);

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

  // Local tensor
  real* A;

  // Local tensor for macro element
  real* macro_A;

  // Array of local dimensions for each argument
  uint* local_dimensions;

  // Array of local dimensions of macro element for primary arguments
  uint* macro_local_dimensions;

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

}

#endif
