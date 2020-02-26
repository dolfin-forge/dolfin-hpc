// Copyright (C) 2006 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-11-29
// Last changed: 2006-11-29

#ifndef __NEW_XML_MESHFUNCTION_H
#define __NEW_XML_MESHFUNCTION_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_XML

#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/io/XMLObject.h>

#include <typeinfo>

namespace dolfin
{

template<typename T>
class XMLMeshFunction : public XMLObject
{
public:

  XMLMeshFunction( MeshFunction<T>& meshfunction );
  ~XMLMeshFunction();

  void startElement (const xmlChar* name, const xmlChar** attrs);
  void endElement   (const xmlChar* name);

  void open(std::string filename);
  bool close();

private:

  enum ParserState { OUTSIDE, INSIDE_MESHFUNCTION, INSIDE_ENTITY, DONE };
  enum MeshFunctionType { INT, UINT, DOUBLE, BOOL, UNSET };

  void readMeshFunction(const xmlChar* name, const xmlChar** attrs);
  void readEntities    (const xmlChar* name, const xmlChar** attrs);

  ParserState state;
  MeshFunctionType mf_type;
  MeshFunction<T> * _meshfunction;

};

//-----------------------------------------------------------------------------
template<typename T>
XMLMeshFunction<T>::XMLMeshFunction(MeshFunction<T>& meshfunction)
  : XMLObject()
  , state(OUTSIDE)
  , mf_type(UNSET)
  , _meshfunction(&meshfunction)
{}
//-----------------------------------------------------------------------------
template<typename T>
XMLMeshFunction<T>::~XMLMeshFunction()
{}
//-----------------------------------------------------------------------------
template<typename T>
void XMLMeshFunction<T>::startElement(const xmlChar* name,
                                      const xmlChar** attrs)
{
  switch ( state )
  {
  case OUTSIDE:
    if ( xmlStrcasecmp(name, (xmlChar *) "meshfunction") == 0 )
    {
      readMeshFunction(name, attrs);
      state = INSIDE_MESHFUNCTION;
    }
    break;
  case INSIDE_MESHFUNCTION:
    if ( xmlStrcasecmp(name, (xmlChar *) "entity") == 0 )
    {
      readEntities(name, attrs);
      state = INSIDE_ENTITY;
    }
    break;
  default:
    break;
  }
}
//-----------------------------------------------------------------------------
template<typename T>
void XMLMeshFunction<T>::endElement(const xmlChar *name)
{
  switch ( state )
  {
  case INSIDE_MESHFUNCTION:
    if ( xmlStrcasecmp(name, (xmlChar *) "meshfunction") == 0 )
    {
      state = DONE;
    }
    break;
  case INSIDE_ENTITY:
    if ( xmlStrcasecmp(name, (xmlChar *) "entity") == 0 )
    {
      state = INSIDE_MESHFUNCTION;
    }
    break;
  default:
    break;
  }
}
//-----------------------------------------------------------------------------
template<typename T>
void XMLMeshFunction<T>::open(std::string filename)
{
  message(1, "Reading mesh from file \"%s\".", filename.c_str());
}
//-----------------------------------------------------------------------------
template<typename T>
bool XMLMeshFunction<T>::close()
{
  return state == DONE;
}
//-----------------------------------------------------------------------------
template<typename T>
void XMLMeshFunction<T>::readMeshFunction(const xmlChar* name,
                                          const xmlChar** attrs)
{
  // Parse values
  std::string type = parse<std::string>(name, attrs, "type");

  if ( type != "int" and type != "uint" and type != "double" and type != "bool" )
    error( "Meshfunction<%s> not implemented", type.c_str() );

  uint tdim = parse<uint>(name, attrs,   "dim");
  uint size = parse<uint>(name, attrs,   "size");

  _meshfunction->init(&(_meshfunction->mesh()), tdim, size);
  if ( _meshfunction->mesh().size(tdim) != size)
    error("Size of mesh function (%d) does not match size of mesh \
          (%d entities of dimension %d)",
          size, _meshfunction->mesh().size(tdim), tdim);

}
//-----------------------------------------------------------------------------
template<typename T>
void XMLMeshFunction<T>::readEntities(const xmlChar* name,
                                      const xmlChar** attrs)
{
  // Read index
  uint i = parse<uint>(name, attrs, "index");
  (*_meshfunction)(i) = parse<T>(name, attrs, "value");
}
//-----------------------------------------------------------------------------

} // namespace dolfin

#endif
#endif
