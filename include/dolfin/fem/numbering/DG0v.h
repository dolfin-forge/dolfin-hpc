// Copyright (C) 2008 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRIVATE_NUMBERING_DG0V_H
#define __DOLFIN_PRIVATE_NUMBERING_DG0V_H

#include <dolfin/fem/DofNumbering.h>

#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

/**
 *  @class  DG0vNumbering
 *
 *  @brief  This class implements a numbering based on the global indices of
 *          cells such that dofs indices are numbered incrementally from offset:
 *
 *            cell index * value dimension
 *
 *          for each component.
 *
 */

class DG0vNumbering : public DofNumbering
{

public:

  ///
  DG0vNumbering(Mesh& mesh, ufc::dofmap& ufc_dofmap) :
      DofNumbering(mesh, ufc_dofmap),
      value_size_(0)
  {
  }

  ///
  ~DG0vNumbering()
  {
  }

  ///
  inline void tabulate_dofs(uint* dofs, ufc::cell const& ufc_cell, Cell const& cell) const
  {
    std::fill_n(dofs, value_size_, value_size_ * ufc_cell.index);
    for (uint k = 1; k < value_size_; ++k)
    {
      ++dofs[k];
    }
  }

  ///
  void build()
  {
    DofNumbering::init();
    //---
    uint const tdim = mesh.topology_dimension();
    Array<ufc::dofmap const*> flattened;
    DofMap::flatten(&ufc_dofmap, flattened);
    value_size_ = flattened.size();
    flattened.free();
    //---
    if (ufc_dofmap.local_dimension() != value_size_)
    {
      error("DG0vNumbering : local dimension %u != %u",
            ufc_dofmap.local_dimension(), value_size_);
    }
    set_range(value_size_ * mesh.topology().offset(tdim),
              value_size_ * mesh.topology().num_owned(tdim));
    //---
    if (mesh.is_distributed())
    {
      DistributedData const& distdata = mesh.distdata()[tdim];
      if (!distdata.valid_numbering)
      {
        error("DG0vNumbering : cell numbering is invalid");
      }
      shared_.clear();
      for (SharedIterator it(distdata); it.valid(); ++it)
      {
        for (uint i = 0; i < value_size_; ++i)
        {
          shared_.insert(value_size_ * it.global_index() + i);
        }
      }
      ghosts_.clear();
      for (GhostIterator it(distdata); it.valid(); ++it)
      {
        for (uint i = 0; i < value_size_; ++i)
        {
          ghosts_.insert(value_size_ * it.global_index() + i);
        }
      }
    }
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
    return std::string("Dof numbering for DG0 vector");
  }

private:

  // Number of scalar entries
  uint value_size_;

  ///
  _set<uint> shared_;
  _set<uint> ghosts_;

};

}
/* namespace dolfin */

#endif /* __DOLFIN_PRIVATE_NUMBERING_DG0V_H */
