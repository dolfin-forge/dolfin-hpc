// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurélien Larcher, 2013.
//
// First added:  2007-04-24
// Last changed: 2013-09-13

#include <dolfin/log/log.h>
#include <dolfin/fem/SubSystem.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
SubSystem::SubSystem()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(uint sub_system)
{
  this->sub_system.push_back(sub_system);
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(uint sub_system, uint sub_sub_system)
{
  this->sub_system.push_back(sub_system);
  this->sub_system.push_back(sub_sub_system);
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(Array<uint> const& sub_system) :
    sub_system(sub_system)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(SubSystem const& sub_system,
                     SubSystem const& sub_sub_system)
{
  this->sub_system = sub_system.sub_system;
  for (uint i = 0; i < sub_sub_system.depth(); ++i)
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
SubSystem const& SubSystem::operator=(SubSystem const& sub_system)
{
  this->sub_system = sub_system.sub_system;
  return *this;
}
//-----------------------------------------------------------------------------
uint SubSystem::depth() const
{
  return sub_system.size();
}
//-----------------------------------------------------------------------------
Array<uint> const& SubSystem::array() const
{
  return sub_system;
}
//-----------------------------------------------------------------------------
std::string SubSystem::str() const
{
  std::stringstream ss;
  ss << "[ ";
  for (uint i = 0; i < sub_system.size(); ++i)
  {
    ss << sub_system[i] << " ";
  }
  ss << "]";
  return ss.str();
}
//-----------------------------------------------------------------------------

}
