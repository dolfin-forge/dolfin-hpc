// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/SubSystem.h>

#include <dolfin/log/log.h>

#include <sstream>

namespace dolfin
{

//-----------------------------------------------------------------------------
SubSystem::SubSystem(size_t sub_system)
{
  this->sub_system.push_back(sub_system);
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(size_t sub_system, size_t sub_sub_system)
{
  this->sub_system.push_back(sub_system);
  this->sub_system.push_back(sub_sub_system);
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(std::vector<size_t> const& sub_system) :
    sub_system(sub_system)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(SubSystem const& sub_system,
                     SubSystem const& sub_sub_system)
{
  this->sub_system = sub_system.sub_system;
  for (size_t i = 0; i < sub_sub_system.depth(); ++i)
  {
    this->sub_system.push_back(sub_sub_system.array()[i]);
  }
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(SubSystem const& sub_system)
{
  this->sub_system = sub_system.sub_system;
}
//-----------------------------------------------------------------------------
auto SubSystem::str() const -> std::string
{
  std::stringstream ss;
  ss << "[ ";
  for ( size_t const & sub : sub_system )
  {
    ss << sub << " ";
  }
  ss << "]";
  return ss.str();
}
//-----------------------------------------------------------------------------

}
