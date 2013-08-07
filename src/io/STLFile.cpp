// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#include <fstream>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/io/STLFile.h>
#include <set>

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
  float data[3];
  uint ntri, v_index, c_index, index[3];
  struct stl_vertex V;
  std::set<stl_vertex> vertices;

  std::ifstream fp(filename.c_str(), std::ifstream::binary);
  fp.read((char *)&hdr, 80*sizeof(char));
  fp.read((char *)&ntri, sizeof(uint));
  
  MeshEditor editor;
  editor.open(mesh, CellType::triangle, 2, 3);  
  editor.initCells(ntri);

  

  v_index = c_index = 0;
  for (uint i = 0; i < ntri; i++) {    
    /* Normal */
    fp.read((char *)&data, 3*sizeof(float)); 

    for (uint j = 0; j < 3; j++) {
      /* Vertex v1 v2 v3 */
      fp.read((char *)&data, 3*sizeof(float));
      V.v1 = (double) data[0];
      V.v2 = (double) data[1];
      V.v3 = (double) data[2];
      
      if (vertices.find(V) == vertices.end()) {
	V.index = v_index++;
	index[j] = V.index;
	vertices.insert(V);

	editor.addVertex(V.index, V.v1, V.v2, V.v3);
      }	
      else {
	index[j] = vertices.find(V)->index;
      }
    }

    editor.addCell(c_index++, index[0], index[1], index[2]);
    
    /* Aux data */
    fp.read((char *)&hdr, 2*sizeof(char));
  }
  fp.close();
  editor.close();

}
//-----------------------------------------------------------------------------


