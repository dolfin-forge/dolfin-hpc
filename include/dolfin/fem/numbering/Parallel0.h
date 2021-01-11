// Copyright (C) 2008 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRIVATE_NUMBERING_PARALLEL0_H
#define __DOLFIN_PRIVATE_NUMBERING_PARALLEL0_H

#include <dolfin/fem/DofNumbering.h>

#include <dolfin/common/timing.h>
#include <dolfin/fem/UFCCellIterator.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/BoundaryMesh.h>

#include <string>

namespace dolfin
{

/**
 *  @class  Parallel0Numbering
 *
 *  @brief  Implements an optimal randomized parallel dof distribution for
 *          scalar elements which does not ensure that all the dofs are owned
 *          by the owner of the entity on which the node is located.
 *
 */

class Parallel0Numbering : public DofNumbering
{

public:

  ///
  Parallel0Numbering(Mesh& mesh, ufc::dofmap& ufc_dofmap) :
      DofNumbering(mesh, ufc_dofmap)
  {
  }

  ///
  ~Parallel0Numbering() override = default;

  ///
  inline void tabulate_dofs(uint* dofs, ufc::cell const&,
                            Cell const& cell) const override
  {
    uint const ii = ufc_dofmap.local_dimension() * cell.index();
    dolfin_assert(array != NULL);
    std::copy(&array[ii], &array[ii] + ufc_dofmap.local_dimension(), dofs);
  }

  ///
  void build() override
  {
    DofNumbering::init();
    //---
    if (!mesh.is_distributed())
    {
      error("Parallel0Numbering : can only be used for a distributed mesh");
    }
    //---

#ifdef HAVE_MPI

    message(1, "Parallel0Numbering : build");
    tic();

    _map<uint, uint> dof_vote;
    _map<uint, std::vector<uint> > dof2index;
    _set<uint> owned;
    _set<uint> ufc_ghosts;
    _set<uint> ufc_shared;
    shared_.clear();
    ghosts_.clear();

    Cell c0(mesh, 0);
    UFCCell ufc_cell(c0);

    uint const pe_size = MPI::size();
    uint const rank = MPI::rank();

    uint const tdim = mesh.topology_dimension();
    uint const local_dim = ufc_dofmap.local_dimension();
    uint * dofs = new uint[local_dim];
    uint const nb_facet_dofs = ufc_dofmap.num_facet_dofs();
    uint * facet_dofs = new uint[nb_facet_dofs];

    // Initialize random number generator differently on each process
    // FIXME: if the ghosts are randomly generated for a given mesh then
    // two instances for the same mesh and same numbering have different
    // ghosts which leads to several issues:
    // - non matching ghosts with consequence of accessing the wrong dof
    // - different local size which leads to crash
    // We have to remove the randomness until a better solution is found.
    // srand((uint)time(0) + MPI::rank());
    std::srand(MPI::seed());

    // List shared dofs and assign vote
    //Array<uint> * sendbuf = new Array<uint>[pe_size];
    std::vector<uint> sendbuf;
    BoundaryMesh& interior_boundary = mesh.interior_boundary();
    for (CellIterator bcell(interior_boundary); !bcell.end(); ++bcell)
    {
      Facet f(mesh, interior_boundary.facet_index(*bcell));
      Cell c(mesh, f.entities(tdim)[0]);

      ufc_cell.update(c);
      ufc_dofmap.tabulate_dofs(dofs, ufc_mesh, ufc_cell);
      uint local_facet = c.index(f);
      ufc_dofmap.tabulate_facet_dofs(facet_dofs, local_facet);

      for (uint i = 0; i < nb_facet_dofs; ++i)
      {
        // Assign an ownership vote for each "shared" dof
        uint dofidx = dofs[facet_dofs[i]];
        if (ufc_shared.find(dofidx) == ufc_shared.end())
        {
          ufc_shared.insert(dofidx);
          uint const vote = (uint) (std::rand() + (real) rank);
          dof_vote[dofidx] = vote;
          sendbuf.push_back(dofidx);
          sendbuf.push_back(vote);
        }
      }
    }

    // Decide ownership of "shared" dofs
    uint src;
    uint dst;
    uint max_recv = sendbuf.size();
    MPI::all_reduce_in_place<MPI::max>( max_recv );
    uint *recvbuf = new uint[max_recv];
    for (uint j = 1; j < pe_size; ++j)
    {
      src = (rank - j + pe_size) % pe_size;
      dst = (rank + j) % pe_size;


      int recv_count = MPI::sendrecv( &sendbuf[0], sendbuf.size(), dst, recvbuf,
                                      max_recv, src, 1 );

      for (int i = 0; i < recv_count; i += 2)
      {
        if (ufc_shared.find(recvbuf[i]) != ufc_shared.end())
        {
          // Move dofs with higher ownership votes from shared to forbidden
          if (recvbuf[i + 1] < dof_vote[recvbuf[i]]
              || (recvbuf[i + 1] == dof_vote[recvbuf[i]] && (src < rank)))
          {
            ufc_ghosts.insert(recvbuf[i]);
            ufc_shared.erase(recvbuf[i]);
          }
        }
      }
    }
    sendbuf.clear();

    // Mark all non forbidden dofs as owned by the processes
    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {
      ufc_cell.update(*cell);
      ufc_dofmap.tabulate_dofs(dofs, ufc_mesh, ufc_cell);
      for (uint i = 0; i < local_dim; ++i)
      {
        if (ufc_ghosts.find(dofs[i]) == ufc_ghosts.end())
        {
          // Mark dof as owned
          owned.insert(dofs[i]);
        }

        // Create mapping from dof to dofmap offset
        dof2index[dofs[i]].push_back(cell->index() * local_dim + i);
      }
    }

    // Cleanup
    delete[] recvbuf;
    delete[] facet_dofs;
    delete[] dofs;

    // Renumber dofs
    array_size = mesh.num_cells() * ufc_dofmap.local_dimension();
    array = new uint[array_size];

    uint const range = owned.size();
    uint local_index = 0;
    MPI::offset(range, local_index);
    uint const rank_offset = local_index;
    set_range(rank_offset, range);

    // Compute renumbering for local and owned shared dofs
    for (_set<uint>::iterator it = owned.begin(); it != owned.end();
         ++it, ++local_index)
    {
      for(std::vector<uint>::iterator di = dof2index[*it].begin();
      di != dof2index[*it].end(); ++di)
      {
        array[*di] = local_index;
      }
      if (ufc_shared.count(*it) > 0)
      {
        sendbuf.push_back(*it);
        sendbuf.push_back(local_index);
        shared_.insert(local_index);
      }
    }
    ufc_shared.clear();
    max_recv = sendbuf.size();
    MPI::all_reduce_in_place<MPI::max>( max_recv );
    recvbuf = new uint[max_recv];
    _set<uint> new_ghosts;
    for (uint j = 1; j < pe_size; ++j)
    {
      src = (rank - j + pe_size) % pe_size;
      dst = (rank + j) % pe_size;

      int recv_count = MPI::sendrecv( &sendbuf[0], sendbuf.size(), dst, recvbuf,
                                      max_recv, src, 1 );

      for (int k = 0; k < recv_count; k += 2)
      {
        // Assign new dof number for ghost dofs
        _set<uint>::iterator it = ufc_ghosts.find(recvbuf[k]);
        if (ufc_ghosts.count(recvbuf[k]) > 0)
        {
          uint const newidx = recvbuf[k + 1];
          dolfin_assert((newidx < rank_offset)||(newidx >= rank_offset + range));
          ufc_ghosts.erase(it);
          // Create ghost set with new numbering
          ghosts_.insert(newidx);
          // Append ghost dofs to owned-shared dofs
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

    tocd();

#endif /* HAVE_MPI */

  }

  ///
  inline bool is_shared(uint index) const override
  {
    return (shared_.count(index) > 0);
  }

  ///
  inline bool is_ghost(uint index) const override
  {
    return (ghosts_.count(index) > 0);
  }

  ///
  inline std::string description() const override
  {
    return std::string("Dof numbering for generic parallel scalar");
  }

private:

  ///
  _set<uint> shared_;
  _set<uint> ghosts_;

};

}
/* namespace dolfin */

#endif /* __DOLFIN_PRIVATE_NUMBERING_PARALLEL0_H */
