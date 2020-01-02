// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/function/SubFunction.h>

#include <dolfin/function/Function.h>
#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Function& SubFunction::function() const
{
  return *f_;
}

//-----------------------------------------------------------------------------
uint SubFunction::index() const
{
  return i_;
}

//-----------------------------------------------------------------------------
void SubFunction::disp() const
{
  section("SubFunction");
  prm("Index", this->index());
  function().disp();
  end();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
