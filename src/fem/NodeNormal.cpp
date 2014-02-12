// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-30
// Last changed: 2014-01-30

#include <dolfin/fem/NodeNormal.h>

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>

#include <map>

namespace dolfin
{

//-----------------------------------------------------------------------------
NodeNormal::NodeNormal(Function& function, VertexNormal::Type weight) :
    mesh_(function.mesh()),
    tdim_(mesh_.topology().dim()),
    normals_(mesh_, weight),
    meshbasis_(normals_.basis()),
    basis_()
{
  if (function.type() == Function::discrete)
  {
    space_ = &(function.space());
    local_space_ = false;
  }
  else
  {
    error("NodeNormal can only be created for discrete functions.");
  }
  for (uint i = 0; i < tdim_; ++i)
  {
    basis_[i].init(mesh_, space_->element().signature());
    V_.push_back(&(basis_[i].vector()));
  }

  ComputeBasisP1();
}

//-----------------------------------------------------------------------------
NodeNormal::~NodeNormal()
{
  Clear();
}

//-----------------------------------------------------------------------------
void NodeNormal::Clear()
{
  V_.clear();
}

//-----------------------------------------------------------------------------
void NodeNormal::ComputeBasisP1()
{
  // Cell tabulated version based on unicorn
  Cell c(mesh_, 0);
  UFCCell ufccell(c);
  uint const local_dim = c.numEntities(0);
  uint *idx = new uint[tdim_ * local_dim];
  uint *id = new uint[tdim_ * local_dim];

  DofMap const& dm = space_->dofmap();
  GenericVector& type = nodetype_.vector();
  MeshDistributedData& distdata = mesh_.distdata();

  // Cell blocks
  real ** cell_block = new real*[tdim_];
  for (uint i = 0; i < tdim_; ++i)
  {
    cell_block[i] = new real[tdim_ * local_dim];
  }
  real * type_block = new real[tdim_ * local_dim];

  // Fill vectors
  for (CellIterator cell(mesh_); !cell.end(); ++cell)
  {

    ufccell.update(*cell, distdata);

    dm.tabulate_dofs(idx, ufccell, cell->index());

    uint ii = 0;
    uint jj = 0;

    // for each space coordinate
    for (uint i = 0; i < tdim_; ++i)
    {
      for (VertexIterator v(*cell); !v.end(); ++v, ++ii)
      {
        if (!distdata.is_ghost(v->index(), 0))
        {
          // for each vector of the basis
          for (uint d = 0; d < tdim_; ++d)
          {
            cell_block[d][jj] = meshbasis_[d]->get(*v);
          }
          type_block[jj] = normals_.vertex_type.get(*v);
          id[jj++] = idx[ii];
        }
      }
    }

    //
    for (uint i = 0; i < tdim_; ++i)
    {
      V_[i]->set(cell_block[i], jj, id);
    }
    type.set(type_block, jj, id);
  }

  for (uint i = 0; i < tdim_; ++i)
  {
    V_[i]->apply();
  }
  type.apply();

  // Cleanup
  delete[] type_block;
  for (uint i = 0; i < tdim_; ++i)
  {
    delete [] cell_block[i];
  }
  delete [] cell_block;
  delete[] id;
  delete[] idx;
}

//-----------------------------------------------------------------------------

}

