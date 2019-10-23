// Copyright (C) 2002-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2002-12-06
// Last changed: 2006-10-16

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_XML

#include <dolfin/io/XMLObject.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
XMLObject::XMLObject()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
XMLObject::~XMLObject()
{
  // Do nothing
}

//-----------------------------------------------------------------------------
void XMLObject::open(std::string const&)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
bool XMLObject::close()
{
  return true;
}
//-----------------------------------------------------------------------------
template <>
const char * XMLObject::strtype(int const&)
{
  return "int";
}
//-----------------------------------------------------------------------------
template <>
int XMLObject::read(const xmlChar * s)
{
  return strtol(reinterpret_cast<const char *>(s), NULL, 0);
}
//-----------------------------------------------------------------------------
template <>
const char * XMLObject::strtype(uint const&)
{
  return "uint";
}
//-----------------------------------------------------------------------------
template <>
uint XMLObject::read(const xmlChar * s)
{
  long value = strtol(reinterpret_cast<const char *>(s), NULL, 0);
  if (value < 0)
  {
    error("XML : parsed unsigned integer is negative");
  }
  return (uint) value;
}
//-----------------------------------------------------------------------------
template <>
const char * XMLObject::strtype(real const&)
{
  return "real";
}
//-----------------------------------------------------------------------------
template <>
real XMLObject::read(const xmlChar * s)
{
  return strtod(reinterpret_cast<const char *>(s), NULL);
}
//-----------------------------------------------------------------------------
template <>
const char * XMLObject::strtype(std::string const&)
{
  return "string";
}
//-----------------------------------------------------------------------------
template <>
std::string XMLObject::read(const xmlChar * s)
{
  return std::string(reinterpret_cast<const char *>(s));
}
//-----------------------------------------------------------------------------
template <>
const char * XMLObject::strtype(bool const&)
{
  return "bool";
}
//-----------------------------------------------------------------------------
template <>
bool XMLObject::read(const xmlChar * s)
{
  const char * value = reinterpret_cast<char const *>(s);
  if (strcmp(value, "true") == 0)
  {
    return true;
  }
  else if (strcmp(value, "false") == 0)
  {
    return false;
  }
  else
  {
    error("Cannot convert \"%s\" to bool.", value);
  }
  return false;
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* HAVE_XML */
