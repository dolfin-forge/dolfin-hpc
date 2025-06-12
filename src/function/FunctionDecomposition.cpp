// Copyright (C) 2010 Niclas Jansson.
// Copyright (C) 2015 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

// Based on Niclas Jansson's AdaptiveRefinement.

#include <dolfin/function/FunctionDecomposition.h>

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/mesh/MeshValues.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>

namespace dolfin
{
//-----------------------------------------------------------------------------
auto FunctionDecomposition::compute(Function const& F) -> std::vector<std::unique_ptr<Function>>
{
  message( 1, "Decomposing Function: %p", &F );

  Mesh& mesh = const_cast<Mesh&>(F.mesh());

  //TODO:
  //This implementation should work on a DiscreteFunction: let us consider
  //that things are acceptable as long as high-level functions are OK.
  FiniteElementSpace const& Wh = F.space();
  ScratchSpace S(Wh);
  std::vector<FiniteElementSpace *> spaces = Wh.flatten();
  std::vector<std::unique_ptr<Function>> Si;
  for ( FiniteElementSpace * fespaces : spaces )
  {
    Si.push_back(std::make_unique<Function>(*fespaces));
  }

  if (Wh.is_vertex_based())
  {
    //NOTE: This implementation is based on the assumption that the dofmap for
    //      a CG1 function is indexed by the global indices of vertices.
    MeshValues<bool, Vertex> marked(mesh);
    size_t offset = 0;
    size_t const numcellnodes = S.local_dimension / S.size;
    real * block = F.create_block();
    F.get_block(block);
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      S.cell.update(*c);
      std::vector<size_t> const & cvi = c->entities(0);
      for (VertexIterator v(*c); !v.end(); ++v)
      {
        size_t ci = 0;
        while (cvi[ci] != v->index())
        {
          ++ci;
        }
        if (!v->is_ghost() && !marked(*v))
        {
          size_t dof_index = v->global_index();
          for (size_t i = 0; i < Si.size(); ++i)
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
    size_t dofii = 0;
    size_t dof_index = 0;
    real * block = F.create_block();
    F.get_block(block);
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      dof_index = c->global_index();
      for ( auto& f : Si )
      {
        f->vector().set(&block[dofii], 1, &dof_index);
        ++dofii;
      }
    }
    delete[] block;
  }
  else if (Wh.element().is_vectorizable())
  {
    size_t ii = 0;
    size_t offset = 0;
    FiniteElementSpace& Vh0 = *(*spaces.begin());
    ScratchSpace S0(Vh0);
    size_t const * celldofs = Vh0.dofmap().dofsmapping();
    real * block = F.create_block();
    F.get_block(block);
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      //
      for ( auto& f : Si )
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
    size_t ii = 0;
    std::vector< size_t > offset( Si.size() );
    std::vector< size_t > local_dim( Si.size() );
    size_t const ** celldofs = new size_t const *[Si.size()];
    for (size_t i = 0; i < Si.size(); ++i)
    {
      FiniteElementSpace const& Vhi = (*Si[i]).space();
      local_dim[i] = Vhi.dofmap().num_element_dofs;
      celldofs[i] = Vhi.dofmap().dofsmapping();
    }
    real * block = F.create_block();
    F.get_block(block);
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      //
      for (size_t i = 0; i < Si.size(); ++i)
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
  for ( auto& f : Si )
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
