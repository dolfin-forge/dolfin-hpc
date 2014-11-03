// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>

#include <fstream>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/io/STLFile.h>
#include <set>

using namespace dolfin;

//-----------------------------------------------------------------------------
STLFile::STLFile(const std::string filename) :
    GenericFile(filename)
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
  uint ntri = 0;
  uint index[3];
  struct stl_vertex V;
  std::set<stl_vertex> vertices;

  std::ifstream fp(filename.c_str(), std::ifstream::binary);
  fp.read((char *) &hdr, 80 * sizeof(char));
  fp.read((char *) &ntri, sizeof(uint));

  MeshEditor editor(mesh, CellType::triangle, 3);
  editor.initCells(ntri);

  uint v_index = 0;
  uint c_index = 0;
  for (uint i = 0; i < ntri; i++)
  {
    /* Normal */
    fp.read((char *) &data, 3 * sizeof(float));

    for (uint j = 0; j < 3; j++)
    {
      /* Vertex v1 v2 v3 */
      fp.read((char *) &data, 3 * sizeof(float));
      V.v[0] = (double) data[0];
      V.v[1] = (double) data[1];
      V.v[2] = (double) data[2];

      if (vertices.find(V) != vertices.end())
      {
        index[j] = vertices.find(V)->index;
      }
      else
      {
        V.index = v_index++;
        index[j] = V.index;
        vertices.insert(V);
      }
    }

    editor.addCell(c_index++, &index[0]);

    /* Aux data */
    fp.read((char *) &hdr, 2 * sizeof(char));
  }

  editor.initVertices(vertices.size());
  for (std::set<stl_vertex>::iterator it = vertices.begin();
      it != vertices.end(); ++it)
  {
    editor.addVertex(it->index, &it->v[0]);
  }

  fp.close();
  editor.close();

}
//-----------------------------------------------------------------------------

