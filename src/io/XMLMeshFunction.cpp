// Copyright (C) 2006 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2007-2008.
//
// First added:  2006-11-27
// Last changed: 2008-04-30

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_XML

#include <cstring>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/io/XMLMeshFunction.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
XMLMeshFunction::XMLMeshFunction(MeshFunction<int>& meshfunction)
: XMLObject(),
  state(OUTSIDE),
  mf_type(UNSET),
  _imeshfunction(&meshfunction),
  _uimeshfunction(0),
  _dmeshfunction(0),
  _bmeshfunction(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
XMLMeshFunction::XMLMeshFunction(MeshFunction<unsigned int>& meshfunction)
: XMLObject(),
  state(OUTSIDE),
  mf_type(UNSET),
  _imeshfunction(0),
  _uimeshfunction(&meshfunction),
  _dmeshfunction(0),
  _bmeshfunction(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
XMLMeshFunction::XMLMeshFunction(MeshFunction<double>& meshfunction)
: XMLObject(),
  state(OUTSIDE),
  mf_type(UNSET),
  _imeshfunction(0),
  _uimeshfunction(0),
  _dmeshfunction(&meshfunction),
  _bmeshfunction(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
XMLMeshFunction::XMLMeshFunction(MeshFunction<bool>& meshfunction)
: XMLObject(),
  state(OUTSIDE),
  mf_type(UNSET),
  _imeshfunction(0),
  _uimeshfunction(0),
  _dmeshfunction(0),
  _bmeshfunction(&meshfunction)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
XMLMeshFunction::~XMLMeshFunction()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void XMLMeshFunction::startElement(const xmlChar *name, const xmlChar **attrs)
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
  default:
    ;
  }
}
//-----------------------------------------------------------------------------
void XMLMeshFunction::endElement(const xmlChar *name)
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
      state = INSIDE_MESHFUNCTION;

    break;

  default:
    ;
  }
}
//-----------------------------------------------------------------------------
void XMLMeshFunction::open(std::string filename)
{
  message(1, "Reading mesh from file \"%s\".", filename.c_str());
}
//-----------------------------------------------------------------------------
bool XMLMeshFunction::close()
{
  return state == DONE;
}
//-----------------------------------------------------------------------------
void XMLMeshFunction::readMeshFunction(const xmlChar *name, const xmlChar **attrs)
{
  // Parse values
  std::string type = parse<std::string>(name, attrs, "type");
  uint tdim = parse<uint>(name, attrs,   "dim");
  uint size = parse<uint>(name, attrs,   "size");
  if (strcmp(type.c_str(), "int") == 0 )
  {
    if (_imeshfunction == NULL)
      error("MeshFunction file of type \"int\", but MeshFunction<int> not initialized");

    _imeshfunction->init(&_uimeshfunction->mesh(), tdim, size);
    if ( _imeshfunction->mesh().size(tdim) != size)
      error("Size of mesh function (%d) does not match size of mesh (%d entities of dimension %d)", size, _imeshfunction->mesh().size(tdim), tdim);
    mf_type = INT;
  }
  else if (strcmp(type.c_str(), "uint") == 0 )
  {
    if (_uimeshfunction == NULL)
      error("MeshFunction file of type \"uint\", but MeshFunction<uint> not initialized");

    _uimeshfunction->init(&_uimeshfunction->mesh(), tdim, size);
    if ( _uimeshfunction->mesh().size(tdim) != size)
      error("Size of mesh function (%d) does not match size of mesh (%d entities of dimension %d)", size, _uimeshfunction->mesh().size(tdim), tdim);
    mf_type = UINT;
  }
  else if (strcmp(type.c_str(), "double") == 0 )
  {
    if (_dmeshfunction == NULL)
      error("MeshFunction file of type \"double\", but MeshFunction<double> not initialized");

    _dmeshfunction->init(&_dmeshfunction->mesh(), tdim, size);
    if ( _dmeshfunction->mesh().size(tdim) != size)
      error("Size of mesh function (%d) does not match size of mesh (%d entities of dimension %d)", size, _dmeshfunction->mesh().size(tdim), tdim);
    mf_type = DOUBLE;
  }
  else if (strcmp(type.c_str(), "bool") == 0 )
  {
    if (_bmeshfunction == NULL)
      error("MeshFunction file of type \"bool\", but MeshFunction<bool> not initialized");

    _bmeshfunction->init(&_uimeshfunction->mesh(), tdim, size);
    if ( _bmeshfunction->mesh().size(tdim) != size)
      error("Size of mesh function (%d) does not match size of mesh (%d entities of dimension %d)", size, _bmeshfunction->mesh().size(tdim), tdim);
    mf_type = BOOL;
  }
}
//-----------------------------------------------------------------------------
void XMLMeshFunction::readEntities(const xmlChar *name, const xmlChar **attrs)
{
  // Read index
  uint i = parse<uint>(name, attrs, "index");
  switch ( mf_type )
  {
  case INT:
    (*_imeshfunction)(i) = parse<int>(name, attrs, "value");
    break;

  case UINT:
    (*_uimeshfunction)(i) = parse<uint>(name, attrs, "value");
    break;

  case DOUBLE:
    (*_dmeshfunction)(i) = parse<real>(name, attrs, "value");
    break;

  case BOOL:
    (*_bmeshfunction)(i) = parse<bool>(name, attrs, "value");
    break;

  default:
    error("Could not determine the type of MeshFunction.");
  }
}
//-----------------------------------------------------------------------------

#endif
