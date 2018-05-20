// Copyright (C) 2008 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __DOLFIN_PRIVATE_NUMBERING_CG1S_H
#define __DOLFIN_PRIVATE_NUMBERING_CG1S_H

#include <dolfin/fem/DofNumbering.h>

#include <dolfin/fem/UFCCell.h>
#include <dolfin/mesh/MeshDistributedData.h>

namespace dolfin
{

/**
 *  @class  CG1sNumbering
 *
 *  @brief  This class implements a numbering based on the global indices of
 *          vertices: tabulating the dofs just amounts to copying the vertex
 *          indices from given UFC cell.
 *          Cell tabulated dofs are not cached.
 *
 */

class CG1sNumbering : public DofNumbering
{

public:

  ///
  CG1sNumbering(Mesh& mesh, ufc::dofmap& ufc_dofmap) :
      DofNumbering(mesh, ufc_dofmap)
  {
  }

  ///
  ~CG1sNumbering()
  {
  }

  ///
  inline void tabulate_dofs(uint* dofs, ufc::cell const& ufc_cell, Cell const& cell) const
  {
    std::copy(ufc_cell.entity_indices[0],
              ufc_cell.entity_indices[0] + cell.num_entities(0), dofs);
  }

  ///
  inline void build()
  {
    DofNumbering::init();
    //---
    if (ufc_dofmap.local_dimension() != mesh.type().num_entities(0))
    {
      error("CG1sNumbering : local dimension %u != %u",
            ufc_dofmap.local_dimension(), mesh.type().num_entities(0));
    }
    set_range(mesh.topology().offset(0), mesh.topology().num_owned(0));
    //---
    if (mesh.is_distributed())
    {
      DistributedData const& distdata = mesh.distdata()[0];
      if (!distdata.valid_numbering)
      {
        error("CG1sNumbering : vertex numbering is invalid");
      }
      shared_.clear();
      for (SharedIterator it(distdata); it.valid(); ++it)
      {
        shared_.insert(it.global_index());
      }
      ghosts_.clear();
      for (GhostIterator it(distdata); it.valid(); ++it)
      {
        ghosts_.insert(it.global_index());
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
    return std::string("Dof numbering for CG1 scalar");
  }

private:

  ///
  _set<uint> shared_;
  _set<uint> ghosts_;

};

}
/* namespace dolfin */

#endif /* __DOLFIN_PRIVATE_NUMBERING_CG1S_H */
