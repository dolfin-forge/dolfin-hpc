// Copyright (C) 2008 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PRIVATE_NUMBERING_DG0S_H
#define __DOLFIN_PRIVATE_NUMBERING_DG0S_H

#include <dolfin/fem/DofNumbering.h>

#include <dolfin/mesh/Mesh.h>

namespace dolfin
{

/**
 *  @class  DG0sNumbering
 *
 *  @brief  This class implements a numbering based on the global indices of
 *          cells. Computation of shared and ghost indices is left but until
 *          support for halo cells is implemented the set are empty.
 *          Cell tabulated dofs are not cached.
 *
 */

class DG0sNumbering : public DofNumbering
{

public:
  ///
  DG0sNumbering( Mesh & mesh, ufc::dofmap & ufc_dofmap )
    : DofNumbering( mesh, ufc_dofmap )
  {
  }

  ///
  ~DG0sNumbering() override = default;

  ///
  inline void tabulate_dofs( size_t *          dofs,
                             ufc::cell const & ufc_cell,
                             Cell const & ) const override
  {
    dofs[0] = ufc_cell.index;
  }

  ///
  void build() override
  {
    DofNumbering::init();
    //---
    size_t const tdim = mesh.topology_dimension();
    if ( ufc_dofmap.num_element_dofs() != 1 )
    {
      error( "DG0sNumbering : local dimension %u != 1",
             ufc_dofmap.num_element_dofs() );
    }
    set_range( mesh.topology().offset( tdim ),
               mesh.topology().num_owned( tdim ) );
    //---
    if ( mesh.is_distributed() )
    {
      DistributedData const & distdata = mesh.distdata()[tdim];
      if ( !distdata.valid_numbering )
      {
        error( "DG0sNumbering : cell numbering is invalid" );
      }
      shared_.clear();
      for ( SharedIterator it( distdata ); it.valid(); ++it )
      {
        shared_.insert( it.global_index() );
      }
      ghosts_.clear();
      for ( GhostIterator it( distdata ); it.valid(); ++it )
      {
        ghosts_.insert( it.global_index() );
      }
    }
  }

  ///
  inline bool is_shared( size_t index ) const override
  {
    return ( shared_.count( index ) > 0 );
  }

  ///
  inline bool is_ghost( size_t index ) const override
  {
    return ( ghosts_.count( index ) > 0 );
  }

  ///
  inline std::string description() const override
  {
    return std::string( "Dof numbering for DG0 scalar" );
  }

private:
  ///
  _set< size_t > shared_;
  _set< size_t > ghosts_;
};

}
/* namespace dolfin */

#endif /* __DOLFIN_PRIVATE_NUMBERING_DG0S_H */
