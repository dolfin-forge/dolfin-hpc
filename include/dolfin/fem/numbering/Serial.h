// Copyright (C) 2015 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRIVATE_NUMBERING_SERIAL_H
#define __DOLFIN_PRIVATE_NUMBERING_SERIAL_H

#include <dolfin/fem/DofNumbering.h>

namespace dolfin
{

/**
 *  @class  SerialNumbering
 *
 */

class SerialNumbering : public DofNumbering
{

public:

  ///
  SerialNumbering(Mesh& mesh, ufc::dofmap& ufc_dofmap) :
      DofNumbering(mesh, ufc_dofmap),
      ufc_mesh_(mesh)
  {
  }

  ///
  ~SerialNumbering() override = default;

  ///
  inline void tabulate_dofs(uint* dofs, ufc::cell const& ufc_cell,
                            Cell const&) const override
  {
    ufc_dofmap.tabulate_dofs(dofs, ufc_mesh_, ufc_cell);
  }

  ///
  inline void build() override
  {
    DofNumbering::init();
    //---
    if (mesh.is_distributed())
    {
      error("SerialNumbering : can only be used for a serial mesh");
    }
    //---
    ufc_mesh_.update();
    set_range(0, ufc_dofmap.global_dimension());
  }

  ///
  inline bool is_shared(uint) const override
  {
    return false;
  }

  ///
  inline bool is_ghost(uint) const override
  {
    return false;
  }

  ///
  inline std::string description() const override
  {
    return std::string("Dof numbering for serial UFC backend");
  }

private:

  ///
  UFCMesh ufc_mesh_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_PRIVATE_NUMBERING_SERIAL_H */
