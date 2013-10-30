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
SubSystem::SubSystem(const Array<uint>& sub_system) :
    sub_system(sub_system)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(const SubSystem& sub_system)
{
  this->sub_system = sub_system.sub_system;
}
//-----------------------------------------------------------------------------
const SubSystem& SubSystem::operator=(const SubSystem& sub_system)
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

}
