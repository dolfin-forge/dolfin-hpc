// Copyright (C) 2006 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-11-29
// Last changed: 2006-11-29

#ifndef __DOLFIN_XML_MESHFUNCTION_H
#define __DOLFIN_XML_MESHFUNCTION_H

#include <dolfin/io/XMLObject.h>

#include <dolfin/mesh/MeshFunction.h>

namespace dolfin
{

template<class T>
class XMLMeshFunction : public XMLObject
{

public:

  ///
  XMLMeshFunction(MeshFunction<T>& meshfunction) :
      XMLObject(),
      state_(ROOT),
      meshfunction_(meshfunction)
  {
    // Do nothing
  }

  ///
  ~XMLMeshFunction()
  {

  }

  ///
  void startElement(const xmlChar* name, const xmlChar** attrs)
  {
    switch (state_)
    {
    case ROOT:
      if (xmlStrcasecmp(name, (xmlChar *) "meshfunction") == 0)
      {
        // Parse type
        std::string type = parse<std::string>(name, attrs, "type");
        T n;
        if (strcmp(type.c_str(), strtype<T>(n)) != 0)
        {
          error("XMLMeshFunction : value type mismatch '%s'", type.c_str());
        }
        // Parse values
        uint tdim = parse<uint>(name, attrs, "dim");
        uint size = parse<uint>(name, attrs, "size");
        uint const meshsize = meshfunction_.mesh().size(tdim);
        if (meshsize != size)
        {
          error("XMLMeshFunction : size mismatch '%u' != '%u'", meshsize, size);
        }
        meshfunction_.init(meshfunction_.mesh(), tdim);
        //
        state_ = IN_MESHFUNCTION;
      }

      break;

    case IN_MESHFUNCTION:

      if (xmlStrcasecmp(name, (xmlChar *) "entity") == 0)
      {
        uint index = parse<uint>(name, attrs, "index");
        meshfunction_.set(index, parse<T>(name, attrs, "value"));
        //
        state_ = IN_ENTITY;
      }
      break;
    default:
      break;
    }
  }

  ///
  void endElement(const xmlChar* name)
  {
    switch (state_)
    {
    case IN_MESHFUNCTION:
      if (xmlStrcasecmp(name, (xmlChar *) "meshfunction") == 0)
      {
        state_ = ROOT;
      }
      break;
    case IN_ENTITY:
      if (xmlStrcasecmp(name, (xmlChar *) "entity") == 0)
      {
        state_ = IN_MESHFUNCTION;
      }
      break;
    default:
      break;
    }
  }

  ///
  void open(std::string const& filename)
  {
    message(1, "Reading mesh from file \"%s\".", filename.c_str());
  }

  ///
  bool close()
  {
    return state_ == ROOT;
  }

private:

  enum ParserState
  {
    ROOT, IN_MESHFUNCTION, IN_ENTITY
  };

  ///
  ParserState state_;
  MeshFunction<T>& meshfunction_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_XML_MESHFUNCTION_H */
