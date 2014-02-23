// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-30
// Last changed: 2014-01-30

#include <dolfin/fem/NodeNormal.h>

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>

#include <map>

namespace dolfin
{

//-----------------------------------------------------------------------------
NodeNormal::NodeNormal(Mesh& mesh, VertexNormal::Type weight) :
    BoundaryNormal(mesh),
    dim_(mesh.geometry().dim()),
    normals_(mesh, weight),
    meshbasis_(normals_.basis())
{
}

//-----------------------------------------------------------------------------
NodeNormal::~NodeNormal()
{
  Clear();
}

//-----------------------------------------------------------------------------
void NodeNormal::compute()
{
  if (V_.size() == 0)
  {
    for (uint i = 0; i < dim_; ++i)
    {
      V_.push_back(&(basis()[i].vector()));
    }
  }

  // Only for Lagrange P1
  ComputeBasisP1();
}

//-----------------------------------------------------------------------------
void NodeNormal::Clear()
{
  V_.clear();
}

//-----------------------------------------------------------------------------
void NodeNormal::ComputeBasisP1()
{
  Mesh& mesh = this->mesh();

  // Cell tabulated version based on unicorn
  Cell c(mesh, 0);
  UFCCell ufccell(c);
  uint const local_dim = c.numEntities(0);
  uint *idx = new uint[dim_ * local_dim];
  uint *id = new uint[dim_ * local_dim];
  uint *type_idx = new uint[local_dim];

  DofMap const& dm = this->basis()[0].space().dofmap();
  DofMap const& dm_type = node_type().space().dofmap();
  GenericVector& type = node_type().vector();
  MeshDistributedData& distdata = mesh.distdata();

  // Cell blocks
  real ** cell_block = new real*[dim_];
  for (uint i = 0; i < dim_; ++i)
  {
    cell_block[i] = new real[dim_ * local_dim];
  }
  real * type_block = new real[local_dim];

  // Fill vectors
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {

    ufccell.update(*cell, distdata);
    dm.tabulate_dofs(idx, ufccell, cell->index());
    dm_type.tabulate_dofs(type_idx, ufccell, cell->index());

    uint ii = 0;
    uint jj = 0;
    // for each space coordinate
    for (uint i = 0; i < dim_; ++i)
    {
      for (VertexIterator v(*cell); !v.end(); ++v, ++ii)
      {
        if (!distdata.is_ghost(v->index(), 0))
        {
          // for each vector of the basis
          for (uint d = 0; d < dim_; ++d)
          {
            cell_block[d][jj] = meshbasis_[d][i].get(*v);
          }
          id[jj++] = idx[ii];
        }
      }
    }
    //
    for (uint d = 0; d < dim_; ++d)
    {
      V_[d]->set(cell_block[d], jj, id);
    }

    //
    for (VertexIterator v(*cell); !v.end(); ++v)
    {
      if (!distdata.is_ghost(v->index(), 0))
      {
        type_block[v.pos()] = normals_.vertex_type.get(*v);
      }
    }
    type.set(type_block, local_dim, type_idx);
  }

  for (uint i = 0; i < dim_; ++i)
  {
    V_[i]->apply();
  }
  type.apply();

  // Cleanup
  delete[] type_block;
  for (uint i = 0; i < dim_; ++i)
  {
    delete[] cell_block[i];
  }
  delete[] cell_block;
  delete[] id;
  delete[] idx;
}

//-----------------------------------------------------------------------------

}

