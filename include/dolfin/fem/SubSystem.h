// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SUB_SYSTEM_H
#define __DOLFIN_SUB_SYSTEM_H

#include <string>
#include <vector>

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
  SubSystem() = default;

  /// Create given sub system (one level)
  SubSystem( size_t sub_system );

  /// Create given sub sub system (two levels)
  SubSystem( size_t sub_system, size_t sub_sub_system );

  /// Create sub system for given array (n levels)
  SubSystem( std::vector< size_t > const & sub_system );

  /// Create given sub and sub sub system (two levels)
  SubSystem( SubSystem const & sub_system, SubSystem const & sub_sub_system );

  /// Copy constructor
  SubSystem( SubSystem const & sub_system );

  /// Assignment operator
  auto operator=( SubSystem const & sub_system ) -> SubSystem & = default;

  /// Return number of levels for nested sub system
  auto depth() const -> size_t;

  /// Return array which defines sub system
  auto array() const -> std::vector< size_t > const &;

  /// Cast
  operator std::vector< size_t > &()
  {
    return sub_system;
  }
  operator std::vector< size_t > const &() const
  {
    return sub_system;
  }

  ///
  auto str() const -> std::string;

private:
  // The array specifying the sub system
  std::vector< size_t > sub_system;
};

//-----------------------------------------------------------------------------
inline auto SubSystem::depth() const -> size_t
{
  return sub_system.size();
}
//-----------------------------------------------------------------------------
inline auto SubSystem::array() const -> std::vector< size_t > const &
{
  return sub_system;
}
//-----------------------------------------------------------------------------

}

#endif
