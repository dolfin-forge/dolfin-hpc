// Copyright (C) 2005-2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/evolution/TimeDependent.h>

#include <dolfin/evolution/Time.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
TimeDependent::TimeDependent(Time const& time) :
    clock_(time.clock()),
    t_(&time)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
TimeDependent& TimeDependent::swap(TimeDependent& other)
{
  std::swap(clock_, other.clock_);
  std::swap(t_, other.t_);
  return *this;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
