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
  RealSpaceNumbering( Mesh & mesh, ufc::dofmap & ufc_dofmap )
    : DofNumbering( mesh, ufc_dofmap )
    , dofs_( NULL )
  {
  }

  ///
  ~RealSpaceNumbering() override
  {
    delete[] dofs_;
  }

  ///
  inline void tabulate_dofs( size_t * dofs,
                             ufc::cell const &,
                             Cell const & ) const override
  {
    std::copy( dofs_, dofs_ + ufc_dofmap.num_element_dofs(), dofs );
  }

  //
  inline void build() override
  {
    DofNumbering::init();

    //---
    size_t offset = ufc_dofmap.num_element_dofs() * MPI::rank();
    set_range( offset, ufc_dofmap.num_element_dofs() );
    delete[] dofs_;
    dofs_ = new size_t[ufc_dofmap.num_element_dofs()];
    std::fill_n( dofs_, ufc_dofmap.num_element_dofs(), offset );
    for ( size_t i = 1; i < ufc_dofmap.num_element_dofs(); ++i )
    {
      dofs_[i] += i;
    }
  }

  ///
  inline bool is_shared( size_t ) const override
  {
    return false;
  }

  ///
  inline bool is_ghost( size_t ) const override
  {
    return false;
  }

  ///
  inline std::string description() const override
  {
    return std::string( "Dof numbering for real space" );
  }

private:
  size_t * dofs_;
};

} /* namespace dolfin */

#endif /* __DOLFIN_PRIVATE_NUMBERING_REALSPACE_H */
