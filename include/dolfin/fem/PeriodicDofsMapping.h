// Copyright (C) 2015 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2015-04-24
// Last changed: 2015-04-24

#ifndef __DOLFIN_PERIODIC_DOFS_MAPPING
#define __DOLFIN_PERIODIC_DOFS_MAPPING

#include <dolfin/common/types.h>

#include <map>

namespace dolfin
{

class DofMap;
class MeshEntity;
class Point;

class PeriodicDofsMapping
{
public:

  ///
  PeriodicDofsMapping(DofMap const& dofmap);

  ///
  ~PeriodicDofsMapping();

  ///
  uint max_local_dimension() const;

  ///
  uint num_Gdofs() const;
  uint num_Hdofs() const;
  uint num_Idofs() const;

  ///
  bool is_Gdof(uint i) const;
  bool is_Hdof(uint i) const;
  bool is_Idof(uint i) const;


  ///
  uint const * get_Gindices() const;

  ///
  void tabulate_dofs(uint Gdof, uint * Hdofs, uint& count) const;

  //
  void tabulate_dofs(uint i, uint * Gdof, uint * Hdofs, uint& count) const;

  ///
  void tabulate_coordinates(uint Gdof, real * Gcoords, real ** Hcoords , uint& count) const;

  //
  void tabulate_coordinates(uint i, uint * Gdof, real * Gcoords, real ** Hcoords, uint& count) const;

  ///
  void disp() const;

private:

  ///
  void init(DofMap const& dofmap);

  ///
  void clear();

  ///
  bool onEntity(Point& p, MeshEntity& entity);

  DofMap const& dofmap_;

  ///
  uint max_local_dimension_;

  /// Map of G dofs to the offset and count in the arrays
  typedef std::map<uint, uint> OffsetMap;
  OffsetMap Goffsets_;
  uint * Gindices_;
  real * Gxcoords_;
  _set<uint> Hdofs_;
  uint * Hcount_;
  uint * Hoffsets_;
  uint * Hindices_;
  real * Hxcoords_;
  _set<uint> Idofs_;

};

}

#endif /* __DOLFIN_PERIODIC_DOFS_MAPPING */
