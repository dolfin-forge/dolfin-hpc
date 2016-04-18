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

class UFCCellIterator : public UFCCell
{

public:

  ///
  UFCCellIterator(Mesh& mesh) :
      UFCCell(),
      it_(mesh),
      geometry_(mesh.geometry())
  {
    cell = &(*it_);
    cell_shape = UFCCell::shape(mesh.type().cellType());
    num_vertices = mesh.type().num_entities(0);
    topological_dimension = mesh.topology().dim();
    geometric_dimension = mesh.geometry().dim();
    // Entities
    entity_indices = new uint*[topological_dimension + 1];
    entity_indices[topological_dimension] = new uint[1];
    for (uint d = 0; d < topological_dimension; ++d)
    {
      entity_indices[d] = new uint[mesh.type().num_entities(d)];
    }
    // Coordinates
    coordinates = new real*[num_vertices];

    //
#if ENABLE_P1_OPTIMIZATIONS
    cell->global_entities(0, entity_indices[0]);
#else
    cell->global_entities(entity_indices);
#endif
    index = entity_indices[topological_dimension][0];
    uint const * vertices = cell->entities(0);
    for (uint i = 0; i < num_vertices; ++i)
    {
      coordinates[i] = geometry_.x(vertices[i]);
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
      it_->global_entities(0, entity_indices[0]);
#else
      it_->global_entities(entity_indices);
#endif
      index = entity_indices[topological_dimension][0];
      uint const * vertices = cell->entities(0);
      for (uint i = 0; i < num_vertices; ++i)
      {
        coordinates[i] = geometry_.x(vertices[i]);
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
  inline Cell* operator->()
  {
    return &(*it_);
  }

  ///
  inline Cell& operator*()
  {
    return (*it_);
  }

private:

  //
  CellIterator it_;

  //
  MeshGeometry& geometry_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_CELL_ITERATOR_H */
