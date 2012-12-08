// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/io/STLFile.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
STLFile::STLFile(const std::string filename) : GenericFile(filename)
{
  type = "STL";
}
//-----------------------------------------------------------------------------
STLFile::~STLFile()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void STLFile::operator>>(Mesh& mesh)
{

}
//-----------------------------------------------------------------------------


