// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-13
// Last changed: 2014-02-13

#ifndef __DOLFIN_SCRATCH_SPACE_H
#define __DOLFIN_SCRATCH_SPACE_H

#include <dolfin/common/types.h>
#include <dolfin/fem/UFCMesh.h>
#include <dolfin/fem/UFCCell.h>

#include <ufc.h>

namespace dolfin
{

class Cell;
class FiniteElementSpace;
class FiniteElement;
class DofMap;
class SubSystem;

/**
 *  @class  ScratchSpace
 *
 *  @brief  Provides data structures for using reference finite element spaces
 *          (i.e on a cell) based on ufc::finite_element and ufc::dofmap.
 */

class ScratchSpace
{
public:

  // Constructor
  ScratchSpace(FiniteElementSpace const& space);

  // Constructor
  ScratchSpace(FiniteElementSpace const& space, SubSystem const& sub_system);

  // Destructor
  ~ScratchSpace();

  // UFC Mesh
  UFCMesh mesh;

  // UFC Cell
  UFCCell cell;

  // Offset for sub system
  uint offset;

  // Finite element
  ufc::finite_element const * finite_element;

  // Dof map
  ufc::dofmap const * dof_map;

  // Value size (number of entries in tensor value)
  uint const size;

  // Reference finite element space dimension
  uint const space_dimension;

  // Reference finite element dof map dimension
  uint const local_dimension;

  // Number of subspaces of the reference finite element
  uint const num_sub_elements;

  // Topological dimension
  uint const topological_dimension;

  // Geometric dimension
  uint const geometric_dimension;

  // Local array for mapping of dofs
  uint * const dofs;

  // Local array for mapping of facet dofs
  uint * const facet_dofs;

  // Local array for values
  real * const values;

  // Local array for expansion coefficients
  real * const coefficients;

  // Local array for basis values
  real * const basis_values;

  // Local array for coordinates
  real** const coordinates;

private:

  uint value_size(ufc::finite_element const& finite_element);

  void init();

  bool const owner_;

};

} /* namespace dolfin */
#endif /* __DOLFIN_SCRATCH_SPACE_H */
