// Copyright (C) 2003-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/common/Variable.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
Variable::Variable()
  : _name( "x" )
  , _label( "data with no label" )
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Variable::Variable( std::string const name, std::string const label )
  : _name( name )
  , _label( label )
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void Variable::rename( std::string const name, std::string const label )
{
  _name  = name;
  _label = label;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
