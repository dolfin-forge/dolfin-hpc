// Copyright (C) 2003-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_VARIABLE_H
#define __DOLFIN_VARIABLE_H

#include <string>

namespace dolfin
{

class Variable
{
public:
  Variable();
  Variable( std::string const name, std::string const label );
  Variable( const Variable & variable ) = default;

  void rename( std::string const name, std::string const label );

  std::string const & name() const;
  std::string const & label() const;

private:
  std::string _name;
  std::string _label;
};

//-----------------------------------------------------------------------------
inline std::string const & Variable::name() const
{
  return _name;
}
//-----------------------------------------------------------------------------
inline std::string const & Variable::label() const
{
  return _label;
}

} // end namespace dolfin

#endif
