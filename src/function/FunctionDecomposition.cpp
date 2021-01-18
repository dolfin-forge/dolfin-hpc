// Copyright (C) 2010 Niclas Jansson.
// Copyright (C) 2015 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

// Based on Niclas Jansson's AdaptiveRefinement.

#include <dolfin/function/FunctionDecomposition.h>

#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/CellIterator.h>
#include <dolfin/mesh/MeshValues.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/VertexIterator.h>

namespace dolfin
{
//-----------------------------------------------------------------------------
Array<Function *> FunctionDecomposition::compute(Function const& F)
{
  message( 1, "Decomposing Function: %p", &F );

  Mesh& mesh = const_cast<Mesh&>(F.mesh());

  //TODO:
  //This implementation should work on a DiscreteFunction: let us consider
  //that things are acceptable as long as high-level functions are OK.
  FiniteElementSpace const& Wh = F.space();
  ScratchSpace S(Wh);
  Array<FiniteElementSpace *> spaces = Wh.flatten();
  Array<Function *> Si;
  for ( FiniteElementSpace * fespaces : spaces )
  {
    Si.push_back(new Function(*fespaces));
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
      Array<uint> const & cvi = c->entities(0);
      for (VertexIterator v(*c); !v.end(); ++v)
      {
        uint ci = 0;
        while (cvi[ci] != v->index())
        {
          ++ci;
        }
        if (!v->is_ghost() && !marked(*v))
        {
          uint dof_index = v->global_index();
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
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      dof_index = c->global_index();
      for ( Function * f : Si )
      {
        f->vector().set(&block[dofii], 1, &dof_index);
        ++dofii;
      }
    }
    delete[] block;
  }
  else if (Wh.element()->is_vectorizable())
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
      for ( Function * f : Si )
      {
        f->vector().set(&block[ii], S0.local_dimension, &celldofs[offset]);
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
    Array< uint > offset( Si.size() );
    Array< uint > local_dim( Si.size() );
    uint const ** celldofs = new uint const *[Si.size()];
    for (uint i = 0; i < Si.size(); ++i)
    {
      FiniteElementSpace const& Vhi = (*Si[i]).space();
      local_dim[i] = Vhi.dofmap().num_element_support_dofs();
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
  }

  // Synchronize leaf functions
  for ( Function * f : Si )
  {
    f->sync();
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
