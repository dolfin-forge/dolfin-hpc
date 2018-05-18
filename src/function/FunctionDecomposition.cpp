// Copyright (C) 2010 Niclas Jansson.
// Copyright (C) 2015 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// Based on Niclas Jansson's AdaptiveRefinement.
//
// First added:
// Last changed:

#include <dolfin/function/FunctionDecomposition.h>

#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/MeshValues.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
FunctionDecomposition::FunctionDecomposition(Function const& F,
                                             Array<Function *>& Si) :
    F_(F),
    Si_(Si)
{
}

//-----------------------------------------------------------------------------
Array<Function *> FunctionDecomposition::compute(Function const& F)
{
  Mesh& mesh = const_cast<Mesh&>(F.mesh());

  //TODO:
  //This implementation should work on a DiscreteFunction: let us consider
  //that things are acceptable as long as high-level functions are OK.
  FiniteElementSpace const& Wh = F.space();
  ScratchSpace S(Wh);
  Array<FiniteElementSpace *> spaces = Wh.flatten();
  Array<Function *> Si;
  for (uint s = 0; s < spaces.size(); ++s)
  {
    Si.push_back(new Function(*(spaces[s])));
  }

  if (Wh.is_vertex_based())
  {
    //NOTE: This implementation is based on the assumption that the dofmap for
    //      a CG1 function is indexed by the global indices of vertices.
    MeshValues<bool, Vertex> marked(mesh);
    uint offset = 0;
    uint const numcellnodes = S.local_dimension / S.size;
    real * block = F.create_block();
    F.get_block(block);
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      S.cell.update(*c);
      uint * cvi = c->entities(0);
      for (VertexIterator v(*c); !v.end(); ++v)
      {
        uint ci = 0;
        while (cvi[ci] != v->index())
        {
          ++ci;
        }
        if (!v->is_ghost() && !marked(*v))
        {
          uint dof_index = mesh.distdata()[0].get_global(v->index());
          for (uint i = 0; i < Si.size(); ++i)
          {
            Si[i]->vector().set(&block[offset + i * numcellnodes + ci], 1,
                                &dof_index);
          }
          marked(*v) = true;
          continue;
        }
      }
      offset += S.local_dimension;
    }
    delete[] block;
  }
  else if (Wh.is_cellwise_constant())
  {
    uint dofii = 0;
    uint dof_index = 0;
    real * block = F.create_block();
    F.get_block(block);
    uint const tdim = mesh.topology_dimension();
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      dof_index = mesh.distdata()[tdim].get_global(c->index());
      for (uint i = 0; i < Si.size(); ++i)
      {
        Si[i]->vector().set(&block[dofii], 1, &dof_index);
        ++dofii;
      }
    }
    delete[] block;
  }
  else if (Wh.element().is_vectorizable())
  {
    uint ii = 0;
    uint offset = 0;
    FiniteElementSpace& Vh0 = *(*spaces.begin());
    ScratchSpace S0(Vh0);
    uint const * celldofs = Vh0.dofmap().dofsmapping();
    real * block = F.create_block();
    F.get_block(block);
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      //
      for (uint i = 0; i < Si.size(); ++i)
      {
        Si[i]->vector().set(&block[ii], S0.local_dimension, &celldofs[offset]);
        ii += S0.local_dimension;
      }
      //
      offset += S0.local_dimension;
    }
    delete[] block;
  }
  else
  {
    uint ii = 0;
    uint * offset = new uint[Si.size()];
    uint * local_dim = new uint[Si.size()];
    uint const ** celldofs = new uint const *[Si.size()];
    for (uint i = 0; i < Si.size(); ++i)
    {
      FiniteElementSpace const& Vhi = (*Si[i]).space();
      local_dim[i] = Vhi.dofmap().local_dimension();
      celldofs[i] = Vhi.dofmap().dofsmapping();
    }
    real * block = F.create_block();
    F.get_block(block);
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      //
      for (uint i = 0; i < Si.size(); ++i)
      {
        Si[i]->vector().set(&block[ii], local_dim[i], &celldofs[i][offset[i]]);
        ii += local_dim[i];
        offset[i] += local_dim[i];
      }
      //
    }
    delete[] block;
    delete[] celldofs;
    delete[] local_dim;
    delete[] offset;
  }

  // Synchronize leaf functions
  for (uint s = 0; s < Si.size(); ++s)
  {
    Si[s]->sync();
  }

  // Cleanup
  while (!spaces.empty())
  {
    delete spaces.back();
    spaces.pop_back();
  }
  spaces.clear();

  return Si;
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
