// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-13
// Last changed: 2014-02-13

#ifndef __SCRATCH_SPACE_H_
#define __SCRATCH_SPACE_H_

#include <dolfin/common/types.h>

#include <dolfin/fem/UFCCell.h>

namespace dolfin
{

class Cell;
class FiniteElement;
class FiniteElementSpace;
class DofMap;

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
  ScratchSpace(Cell& cell, FiniteElement const& finite_element,
               DofMap const& dof_map);

  // Constructor
  ScratchSpace(FiniteElementSpace const& space);

  // Destructor
  ~ScratchSpace();

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

  // Tabulation on reference cell
  uint * const & cell_tabulation() const;

  // Tabulation per subspace
  uint ** const & sub_element_cell_tabulation() const;

  // Tabulation per entity
  uint ** const & entity_cell_tabulation() const;

  // Local array for mapping of dofs
  uint * const dofs;

  // Local array for expansion coefficients
  real * const coefficients;

  // Local array for values
  real * const values;

  // Local array for coordinates
  real** const coordinates;

  // UFC Cell
  UFCCell cell;

private:

  void Initialize(Cell const& cell, FiniteElement const& finite_element,
                  DofMap const& dof_map);

  uint value_size(FiniteElement const& finite_element);

  void set_cell_tabulation(Cell const& cell, DofMap const& dof_map,
                           uint **& dofs);

  uint * const tabulation_on_cell_;

  uint ** const tabulation_per_sub_element_;

  uint * const sub_element_space_dimensions_;

  uint ** const tabulation_per_entity_;

  uint * const num_entity_dofs_;

};

} /* namespace dolfin */
#endif /* __SCRATCH_SPACE_H_ */
