// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#include <fstream>
#include <dolfin/mesh/MeshEditor.h>
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

  char hdr[80];
  uint ntri;


  std::ifstream fp(filename.c_str(), std::ifstream::binary);
  fp.read((char *)&hdr, 80*sizeof(char));
  fp.read((char *)&ntri, sizeof(uint));
  
  MeshEditor editor;
  editor.open(mesh, CellType::triangle, 2, 3);  
  editor.close();

  for (uint i = 0; i < ntri; i++) {
    /*
      REAL32[3] Normal vector
      REAL32[3]  Vertex 1
      REAL32[3]  Vertex 2
      REAL32[3]  Vertex 3
      UINT16  Attribute byte count
    */
  }


  fp.close();
}
//-----------------------------------------------------------------------------


