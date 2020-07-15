// Copyright (C) 2008 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRIVATE_NUMBERING_CG1V_H
#define __DOLFIN_PRIVATE_NUMBERING_CG1V_H

#include <dolfin/fem/DofNumbering.h>

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/mesh/MeshDistributedData.h>
#include <dolfin/mesh/Vertex.h>

#include <string>

namespace dolfin
{

/**
 *  @class  CG1vNumbering
 *
 *  @brief  This class implements a numbering based on the global indices of
 *          vertices. The implementation exchanging degrees of freedom to
 *          renumber owned degrees of freedom contiguously per process is
 *          removed for now since vertices are always renumbered that way.
 *          Cell tabulated dofs are not cached.
 *
 */

class CG1vNumbering : public DofNumbering
{

public:

  ///
  CG1vNumbering(Mesh& mesh, ufc::dofmap& ufc_dofmap) :
      DofNumbering(mesh, ufc_dofmap),
      value_size_(0)
  {
  }

  ///
  ~CG1vNumbering()
  {
  }

  /// Tabulate dofs on cell
  inline void tabulate_dofs(uint* dofs, ufc::cell const& ufc_cell, Cell const& cell) const
  {
    for (uint k = 0; k < value_size_; ++k)
    {
      std::copy(ufc_cell.entity_indices[0],
                ufc_cell.entity_indices[0] + cell.num_entities(0), dofs);
      for (uint v = 0; v < cell.num_entities(0); ++v, ++dofs)
      {
        *dofs *= value_size_;
        *dofs += k;
      }
    }
  }

  /// Build dofmap
  void build()
  {
    DofNumbering::init();
    //---
    Array<ufc::dofmap const*> flattened;
    DofMap::flatten(&ufc_dofmap, flattened);
    value_size_ = flattened.size();
    free( flattened );
    //---
    if (ufc_dofmap.local_dimension()
        != mesh.type().num_entities(0) * value_size_)
    {
      error("CG1sNumbering : local dimension %u != %u",
            ufc_dofmap.local_dimension(),
            mesh.type().num_entities(0) * value_size_);
    }
    set_range(value_size_ * mesh.topology().offset(0),
              value_size_ * mesh.topology().num_owned(0));
    //---
    if (mesh.is_distributed())
    {
      DistributedData const& distdata = mesh.distdata()[0];
      if (!distdata.valid_numbering)
      {
        error("CG1vNumbering : vertex numbering is invalid");
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
    return std::string("Dof numbering for CG1 vector");
  }

private:

  /// Number of scalar entries
  uint value_size_;

  ///
  _set<uint> shared_;
  _set<uint> ghosts_;

};

}
/* namespace dolfin */

#endif /* __DOLFIN_PRIVATE_NUMBERING_CG1V_H */
