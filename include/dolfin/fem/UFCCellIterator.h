// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __DOLFIN_UFC_CELL_ITERATOR_H
#define __DOLFIN_UFC_CELL_ITERATOR_H

#include <dolfin/common/types.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/mesh/Cell.h>

namespace dolfin
{

class UFCCellIterator
{

public:

  ///
  UFCCellIterator(Mesh& mesh) :
      ufc_cell_(),
      it_(mesh),
      geometry_(mesh.geometry())
  {
    ufc_cell_.cell_ = &(*it_);
    ufc_cell_.cell_shape = UFCCell::shape(mesh.type().cellType());
    ufc_cell_.num_vertices = mesh.type().num_entities(0);
    ufc_cell_.topological_dimension = mesh.topology().dim();
    ufc_cell_.geometric_dimension = mesh.geometry().dim();
    // Entities
    ufc_cell_.entity_indices = new uint*[ufc_cell_.topological_dimension + 1];
    ufc_cell_.entity_indices[ufc_cell_.topological_dimension] = new uint[1];
    for (uint d = 0; d < ufc_cell_.topological_dimension; ++d)
    {
      ufc_cell_.entity_indices[d] = new uint[mesh.type().num_entities(d)];
    }
    // Coordinates
    ufc_cell_.coordinates = new real*[ufc_cell_.num_vertices];

    //
#if ENABLE_P1_OPTIMIZATIONS
    ufc_cell_.cell_->get_global_entities(0, ufc_cell_.entity_indices[0]);
#else
    ufc_cell_.cell_->get_global_entities(ufc_cell_.entity_indices);
#endif
    ufc_cell_.index =
        ufc_cell_.entity_indices[ufc_cell_.topological_dimension][0];
    uint const * vertices = ufc_cell_.cell_->entities(0);
    for (uint i = 0; i < ufc_cell_.num_vertices; ++i)
    {
      ufc_cell_.coordinates[i] = geometry_.x(vertices[i]);
    }
  }

  ///
  ~UFCCellIterator()
  {
  }

  inline UFCCellIterator& operator++()
  {
    if (!(++it_).end())
    {
      // The mesh iterator should be used instead of the cell pointer because
      // the underlying mesh entity is only updated when dereferenced from
      // iterator.
#if ENABLE_P1_OPTIMIZATIONS
      it_->get_global_entities(0, ufc_cell_.entity_indices[0]);
#else
      it_->get_global_entities(ufc_cell_.entity_indices);
#endif
      ufc_cell_.index =
          ufc_cell_.entity_indices[ufc_cell_.topological_dimension][0];
      uint const * vertices = ufc_cell_.cell_->entities(0);
      for (uint i = 0; i < ufc_cell_.num_vertices; ++i)
      {
        ufc_cell_.coordinates[i] = geometry_.x(vertices[i]);
      }
    }
    return *this;
  }

  ///
  inline bool end() const
  {
    return it_.end();
  }

  ///
  inline UFCCell* operator->()
  {
    return &(ufc_cell_);
  }

  ///
  inline UFCCell& operator*()
  {
    return (ufc_cell_);
  }

  ///
  inline Cell& cell()
  {
    return (*it_);
  }

private:

  //
  UFCCell ufc_cell_;
  CellIterator it_;

  //
  MeshGeometry& geometry_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_CELL_ITERATOR_H */
