// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SCRATCH_SPACE_H
#define __DOLFIN_SCRATCH_SPACE_H

#include <dolfin/common/types.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/ufc/ufc.h>

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

  // UFC Cell
  UFCCell cell;

  // Offset for sub system
  size_t offset;

  // Finite element
  ufc::finite_element const * finite_element;

  // Dof map
  ufc::dofmap const * dof_map;

  // Value size (number of entries in tensor value)
  size_t const size;

  // Reference finite element space dimension
  size_t const space_dimension;

  // Reference finite element dof map dimension
  size_t const local_dimension;

  // Number of subspaces of the reference finite element
  size_t const num_sub_elements;

  // Topological dimension
  size_t const topological_dimension;

  // Geometric dimension
  size_t const geometric_dimension;

  // Local array for mapping of dofs
  size_t * const dofs;

  // Local array for mapping of facet dofs
  size_t * const facet_dofs;

  // Local array for values
  real * const values;

  // Local array for expansion coefficients
  real * const coefficients;

  // Local array for basis values
  real * const basis_values;

#ifdef ENABLE_EVALUATE_BASIS_FROM_COORDINATES
  // Local array for all basis values
  real ** const all_basis_values;
#endif

  // Local array for coordinates
  std::vector< real > coordinates;

private:

  // Copy constructor
  ScratchSpace(ScratchSpace const& other);

  /// FIXME this should now be present in the ufc::FE interface
  auto value_size(ufc::finite_element const& finite_element) -> size_t;

  void init();

  bool const owner_;

};

} /* namespace dolfin */
#endif /* __DOLFIN_SCRATCH_SPACE_H */
