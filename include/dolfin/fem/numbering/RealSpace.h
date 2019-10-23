// Copyright (C) 2014 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRIVATE_NUMBERING_REALSPACE_H
#define __DOLFIN_PRIVATE_NUMBERING_REALSPACE_H

#include <dolfin/fem/DofNumbering.h>

namespace dolfin
{

class RealSpaceNumbering : public DofNumbering
{

public:

  ///
  RealSpaceNumbering(Mesh& mesh, ufc::dofmap& ufc_dofmap) :
      DofNumbering(mesh, ufc_dofmap),
      dofs_(NULL)
  {
  }

  ///
  ~RealSpaceNumbering()
  {
    delete [] dofs_;
  }

  ///
  inline void tabulate_dofs(uint* dofs, ufc::cell const& ufc_cell, Cell const& cell) const
  {
    std::copy(dofs_, dofs_ + ufc_dofmap.local_dimension(), dofs);
  }

  //
  inline void build()
  {
    DofNumbering::init();

    //---
    uint offset = ufc_dofmap.local_dimension() * MPI::rank();
    set_range(offset, ufc_dofmap.local_dimension());
    delete [] dofs_;
    dofs_ = new uint[ufc_dofmap.local_dimension()];
    std::fill_n(dofs_, ufc_dofmap.local_dimension(), offset);
    for (uint i = 1; i < ufc_dofmap.local_dimension(); ++i)
    {
      dofs_[i] += i;
    }
  }

  ///
  inline bool is_shared(uint index) const
  {
    return false;
  }

  ///
  inline bool is_ghost(uint index) const
  {
    return false;
  }

  ///
  inline std::string description() const
  {
    return std::string("Dof numbering for real space");
  }

private:

  uint * dofs_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_PRIVATE_NUMBERING_REALSPACE_H */
