// Copyright (C) 2013 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRIVATE_NUMBERING_PARALLEL1_H
#define __DOLFIN_PRIVATE_NUMBERING_PARALLEL1_H

#include <dolfin/fem/DofNumbering.h>
#include <dolfin/common/DistributedData.h>

#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/BoundaryMesh.h>

#include <cstring>

namespace dolfin
{

/**
 *  @class  Parallel1Numbering
 *
 *  @brief  Implements a parallel dof distribution for vector elements which
 *          ensures that the ownership of dofs matches the ownership of mesh
 *          entities on which they are located.
 *
 */

class Parallel1Numbering : public DofNumbering
{

public:

  ///
  Parallel1Numbering(Mesh& mesh, ufc::dofmap& ufc_dofmap) :
      DofNumbering(mesh, ufc_dofmap)
  {
  }

  ///
  ~Parallel1Numbering()
  {
  }

  ///
  inline void tabulate_dofs(uint* dofs, ufc::cell const&,
                            Cell const& cell) const
  {
    uint const ii = ufc_dofmap.local_dimension() * cell.index();
    std::copy(&array[ii], &array[ii] + ufc_dofmap.local_dimension(), dofs);
  }

  ///
  void build()
  {
    DofNumbering::init();
    //---
    if (!mesh.is_distributed())
    {
      error("Parallel1Numbering : can only be used for a distributed mesh");
    }
    //---

#ifdef HAVE_MPI

    /*
     * This is a naive and inefficient renumbering hastly coded, sorry...
     */

    UFCMesh ufc_mesh(mesh);
    uint const tdim = mesh.topology_dimension();
    uint * dofs = new uint[ufc_dofmap.local_dimension()];
    _set<uint> owned;
    _set<uint> ufc_ghosts;
    _set<uint> ufc_shared;
    _map<uint, std::vector<uint> > dof2index;
    shared_.clear();
    ghosts_.clear();

    uint pe_size = MPI::size();
    uint rank = MPI::rank();

    uint num_expected_shared = 0;
    uint num_expected_ghosts = 0;
    for (uint i = 0; i < tdim; ++i)
    {
      num_expected_shared += ufc_dofmap.num_entity_dofs(i)
          * mesh.topology().num_shared(i);
      num_expected_ghosts += ufc_dofmap.num_entity_dofs(i)
          * mesh.topology().num_ghost(i);
    }

    // Cache number of dofs per mesh entity
    uint * num_entity_dofs = new uint[tdim + 1];
    uint max_num_entity_dof = 0;
    for (uint d = 0; d <= tdim; ++d)
    {
      num_entity_dofs[d] = ufc_dofmap.num_entity_dofs(d);
      max_num_entity_dof = std::max(max_num_entity_dof, num_entity_dofs[d]);
    }
    dolfin_assert(max_num_entity_dof > 0);
    uint * entity_local_dofs = new uint[max_num_entity_dof];
    uint * entity_dofs = new uint[max_num_entity_dof];
    uint * facet_dofs = new uint[ufc_dofmap.num_facet_dofs()];

    // Loop, baby, loop !
    uint ii = 0;
    CellIterator cell(mesh);
    UFCCell ufc_cell(*cell);
    for (; !cell.end(); ++cell)
    {
      ufc_cell.update(*cell);
      ufc_dofmap.tabulate_dofs(dofs, ufc_mesh, ufc_cell);

      // Create mapping from dof to dofmap offset
      for (uint i = 0; i < ufc_dofmap.local_dimension(); ++i)
      {
        dof2index[dofs[i]].push_back(ii++);
      }

      // Add dofs restricted to the cell as owned
      uint const num_celldofs = num_entity_dofs[tdim];
      ufc_dofmap.tabulate_entity_dofs(entity_local_dofs, tdim, 0);
      for (uint dof = 0; dof < num_celldofs; ++dof)
      {
        owned.insert(dofs[entity_local_dofs[dof]]);
      }

      // Decide ownership
      for (uint d = 0; d < tdim; ++d)
      {
        uint const num_dofs = num_entity_dofs[d];
        for (MeshEntityIterator m(*cell, d); !m.end(); ++m)
        {
          dolfin_assert(m.pos() == (uint)cell->index(*m));
          // Get the dof indices for the entity
          ufc_dofmap.tabulate_entity_dofs(entity_local_dofs, d, m.pos());
          for (uint dof = 0; dof < num_dofs; ++dof)
          {
            entity_dofs[dof] = dofs[entity_local_dofs[dof]];
          }

          //
          if (m->is_ghost())
          {
            dolfin_assert(m->is_shared());
            // Append UFC indices of ghost entities
            ufc_ghosts.insert(entity_dofs, entity_dofs + num_dofs);
          }
          else if (m->is_shared())
          {
            // Append owned shared indices: these will not be renumbered
            ufc_shared.insert(entity_dofs, entity_dofs + num_dofs);
            owned.insert(entity_dofs, entity_dofs + num_dofs);
          }
          else
          {
            owned.insert(entity_dofs, entity_dofs + num_dofs);
          }
        }
      }
    }
    delete[] facet_dofs;
    delete[] entity_dofs;
    delete[] entity_local_dofs;
    delete[] num_entity_dofs;
    delete[] dofs;

    dolfin_assert(ufc_ghosts.size() == num_expected_ghosts);
    dolfin_assert(ufc_shared.size() == num_expected_shared - num_expected_ghosts);

    // Renumber dofs
    array_size = mesh.num_cells() * ufc_dofmap.local_dimension();
    array = new uint[array_size];

    uint const range = owned.size();
    uint local_index = 0;
    MPI::offset(range, local_index);
    uint const rank_offset = local_index;
    set_range(rank_offset, range);

    // Compute renumbering for local and owned shared dofs
    uint src, dst, max_recv;
    std::vector<uint> sendbuf;
    for (_set<uint>::iterator it = owned.begin(); it != owned.end();
         ++it, ++local_index)
    {
      for(std::vector<uint>::iterator di = dof2index[*it].begin();
      di != dof2index[*it].end(); ++di)
      {
        array[*di] = local_index;
      }
      // Send owned shared dofs: old index then new index
      if (ufc_shared.count(*it) > 0)
      {
        sendbuf.push_back(*it);
        sendbuf.push_back(local_index);
        shared_.insert(local_index);
      }
    }
    ufc_shared.clear();
    uint local_size = sendbuf.size();
    MPI::all_reduce<MPI::max>( local_size, max_recv );
    uint * recvbuf = new uint[max_recv];
    for (uint j = 1; j < pe_size; ++j)
    {
      src = (rank - j + pe_size) % pe_size;
      dst = (rank + j) % pe_size;

      int recv_count = MPI::sendrecv( &sendbuf[0], sendbuf.size(), dst,
                                      recvbuf, max_recv, src,1 );

      for (int k = 0; k < recv_count; k += 2)
      {
        // Assign new dof number for ghost dofs
        _set<uint>::iterator it = ufc_ghosts.find(recvbuf[k]);
        if (it != ufc_ghosts.end())
        {
          uint const newidx = recvbuf[k + 1];
          dolfin_assert((newidx < rank_offset)||(newidx >= rank_offset + range));
          ufc_ghosts.erase(it);
          ghosts_.insert(newidx);
          dolfin_assert(shared_.count(newidx) == 0);
          shared_.insert(newidx);
          for (std::vector<uint>::iterator di = dof2index[recvbuf[k]].begin();
               di != dof2index[recvbuf[k]].end(); ++di)
          {
            array[*di] = newidx;
          }
        }
      }
    }
    delete[] recvbuf;

    message(1, "Parallel1Numbering: Remaining ghost dofs %u", ufc_ghosts.size());

    //---

    if (ghosts_.size() != num_expected_ghosts)
    {
      error("Mismatch: expected number of shared dofs and actual, %u != %u",
            num_expected_ghosts, ghosts_.size());
    }
    if (shared_.size() != num_expected_shared)
    {
      error("Mismatch: expected number of shared dofs and actual, %u != %u",
            num_expected_shared, shared_.size());
    }

#endif /* HAVE_MPI */

  }

  ///
  inline bool is_shared(uint index) const
  {
    return (shared_.count(index) > 0);
  }

  ///
  inline bool is_ghost(uint index) const
  {
    return (ghosts_.count(index) > 0);
  }

  ///
  inline std::string description() const
  {
    return std::string("Dof numbering for generic parallel vector");
  }

private:

  ///
  _set<uint> shared_;
  _set<uint> ghosts_;

};

}
/* namespace dolfin */

#endif /* __DOLFIN_PRIVATE_NUMBERING_PARALLEL1_H */
