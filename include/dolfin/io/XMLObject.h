// Copyright (C) 2003-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_XML_OBJECT_H
#define __DOLFIN_XML_OBJECT_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_XML

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

#include <libxml/parser.h>

#include <cstring>
#include <string>

namespace dolfin
{

class XMLObject
{

public:

  /// Constructor
  XMLObject();

  /// Destructor
  virtual ~XMLObject();

  /// Callback for start of XML element
  virtual void startElement(const xmlChar* name, const xmlChar** attrs) = 0;

  /// Callback for end of XML element
  virtual void endElement(const xmlChar* name) = 0;

  /// Callback for start of XML file (optional)
  virtual void open(std::string const& filename);

  /// Callback for end of XML file, should return true iff data is ok (optional)
  virtual bool close();

protected:

  // Template function for value type
  template <class T>
  const char * strtype(T const&)
  {
    error("XMLObject : string type unimplemented for requested value type");
    return T();
  }

  // Template function for reading values
  template <class T>
  T read(const xmlChar *)
  {
    error("XMLObject : reading function unimplemented for given type");
    return T();
  }

  // Template function for parsing
  template <class T>
  T parse(const xmlChar* name, const xmlChar** attrs, const char * label)
  {
    if (attrs == NULL)
    {
      error("Missing attribute \"%s\" for <%s> in XML file.", label, name);
    }
    for (uint i = 0; attrs[i] != NULL; ++i)
    {
      if (xmlStrcasecmp(attrs[i], (xmlChar *) label) == 0)
      {
        if (attrs[i + 1] == NULL)
        {
          error("Value for attribute \"%s\" of <%s> missing in XML file.",
                label, name);
        }
        return read<T>(attrs[i + 1]);
      }
    }
    error("Missing attribute \"%s\" for <%s> in XML file.", label, name);
    return T();
  }

};

//--- TEMPLATES ---------------------------------------------------------------

template <>
const char * XMLObject::strtype(int const& t);

//-----------------------------------------------------------------------------
template <>
int XMLObject::read(const xmlChar * s);

//-----------------------------------------------------------------------------
template <>
const char * XMLObject::strtype(uint const& t);

//-----------------------------------------------------------------------------
template <>
uint XMLObject::read(const xmlChar * s);

//-----------------------------------------------------------------------------
template <>
const char * XMLObject::strtype(real const& t);

//-----------------------------------------------------------------------------
template <>
real XMLObject::read(const xmlChar * s);

//-----------------------------------------------------------------------------
template <>
const char * XMLObject::strtype(std::string const& t);

//-----------------------------------------------------------------------------
template <>
std::string XMLObject::read(const xmlChar * s);

//-----------------------------------------------------------------------------
template <>
const char * XMLObject::strtype(bool const& t);

//-----------------------------------------------------------------------------
template <>
bool XMLObject::read(const xmlChar * s);

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* HAVE_XML */

#endif /* __DOLFIN_XML_OBJECT_H */
