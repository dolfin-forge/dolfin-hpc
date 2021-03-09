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

  auto name() const -> std::string const &;
  auto label() const -> std::string const &;

private:
  std::string _name;
  std::string _label;
};

//-----------------------------------------------------------------------------
inline auto Variable::name() const -> std::string const &
{
  return _name;
}
//-----------------------------------------------------------------------------
inline auto Variable::label() const -> std::string const &
{
  return _label;
}

} // end namespace dolfin

#endif
