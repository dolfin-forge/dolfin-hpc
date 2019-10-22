// Copyright (C) 2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2005-12-18
// Last changed: 2005-12-21

#include <dolfin/log/dolfin_log.h>
#include <dolfin/parameter/ParameterValue.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
ParameterValue::ParameterValue()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ParameterValue::~ParameterValue()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
ParameterValue const& ParameterValue::operator=(int)
{
  error("Cannot assign int value to parameter of type %s.", type().c_str());
  return *this;
}
//-----------------------------------------------------------------------------
ParameterValue const& ParameterValue::operator=(real)
{
  error("Cannot assign real value to parameter of type %s.", type().c_str());
  return *this;
}
//-----------------------------------------------------------------------------
ParameterValue const& ParameterValue::operator=(bool)
{
  error("Cannot assign bool value to parameter of type %s.", type().c_str());
  return *this;
}
//-----------------------------------------------------------------------------
ParameterValue const& ParameterValue::operator=(std::string)
{
  error("Cannot assign string value to parameter of type %s.", type().c_str());
  return *this;
}
//-----------------------------------------------------------------------------
ParameterValue const& ParameterValue::operator=(uint)
{
  error("Cannot assign uint value to parameter of type %s.", type().c_str());
  return *this;
}
//-----------------------------------------------------------------------------
ParameterValue::operator int() const
{
  error("Unable to convert parameter of type %s to int.", type().c_str());
  return 0;
}
//-----------------------------------------------------------------------------
ParameterValue::operator real() const
{
  error("Unable to convert parameter of type %s to real.", type().c_str());
  return 0.0;
}
//-----------------------------------------------------------------------------
ParameterValue::operator bool() const
{
  error("Unable to convert parameter of type %s to bool.", type().c_str());
  return false;
}
//-----------------------------------------------------------------------------
ParameterValue::operator std::string() const
{
  error("Unable to convert parameter of type %s to string.", type().c_str());
  return "";
}
//-----------------------------------------------------------------------------
ParameterValue::operator uint() const
{
  error("Unable to convert parameter of type %s to uint.", type().c_str());
  return 0;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
