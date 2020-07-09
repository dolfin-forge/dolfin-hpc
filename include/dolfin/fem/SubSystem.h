// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SUB_SYSTEM_H
#define __DOLFIN_SUB_SYSTEM_H

#include <dolfin/common/Array.h>

namespace dolfin
{

class Mesh;

/// This class represents a sub system that may be specified as a
/// recursively nested sub system of some given system.
///
/// The sub system is specified by an array of indices. For example,
/// the array [3, 0, 2] specifies sub system 2 of sub system 0 of
/// sub system 3.

class SubSystem
{
public:
  /// Create empty sub system (no sub systems)
  SubSystem();

  /// Create given sub system (one level)
  SubSystem( uint sub_system );

  /// Create given sub sub system (two levels)
  SubSystem( uint sub_system, uint sub_sub_system );

  /// Create sub system for given array (n levels)
  SubSystem( Array< uint > const & sub_system );

  /// Create given sub and sub sub system (two levels)
  SubSystem( SubSystem const & sub_system, SubSystem const & sub_sub_system );

  /// Copy constructor
  SubSystem( SubSystem const & sub_system );

  /// Assignment operator
  SubSystem const & operator=( SubSystem const & sub_system );

  /// Return number of levels for nested sub system
  uint depth() const;

  /// Return array which defines sub system
  Array< uint > const & array() const;

  /// Cast
  operator Array< uint > &()
  {
    return sub_system;
  }
  operator Array< uint > const &() const
  {
    return sub_system;
  }

  ///
  std::string str() const;

private:
  // The array specifying the sub system
  Array< uint > sub_system;
};

//-----------------------------------------------------------------------------
inline SubSystem const & SubSystem::operator=( SubSystem const & sub_system )
{
  this->sub_system = sub_system.sub_system;
  return *this;
}
//-----------------------------------------------------------------------------
inline uint SubSystem::depth() const
{
  return sub_system.size();
}
//-----------------------------------------------------------------------------
inline Array< uint > const & SubSystem::array() const
{
  return sub_system;
}
//-----------------------------------------------------------------------------

}

#endif
