// Copyright (C) 2012 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/io/STLFile.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/byteswap.h>
#include <dolfin/mesh/MeshEditor.h>

#include <fstream>

namespace dolfin
{

//--- Implementation detail ---------------------------------------------------
struct stl_vertex
{
  double v[3];
  uint index;

  bool operator <(const stl_vertex& other) const
  {
    return ((v[0] < other.v[0])
        || (v[0] == other.v[0]
            && (v[1] < other.v[1] || (v[1] == other.v[1] && v[2] < other.v[2]))));
  }

  bool operator ==(const stl_vertex& other) const
  {
    return (v[0] == other.v[0] && v[1] == other.v[1] && v[2] == other.v[2]);
  }
};

//-----------------------------------------------------------------------------
STLFile::STLFile(const std::string filename) :
    GenericFile("STL", filename)
{
}
//-----------------------------------------------------------------------------
void STLFile::operator>>(Mesh& mesh)
{
  char hdr[80];
  float data[3];
  uint ntri = 0;
  uint index[3];
  struct stl_vertex V;
  _ordered_set<stl_vertex> vertices;

  std::ifstream fp(filename.c_str(), std::ifstream::binary);
  fp.read((char *) &hdr, 80 * sizeof(char));
  fp.read((char *) &ntri, sizeof(uint));

#ifdef HAVE_BIG_ENDIAN
  ntri = bswap(ntri);
#endif

  MeshEditor editor(mesh, CellType::triangle, 3, dolfin::MPI::DOLFIN_COMM_SELF);
  editor.init_cells(ntri);

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

#ifdef HAVE_BIG_ENDIAN
      V.v[0] = (double) bswap(data[0]);
      V.v[1] = (double) bswap(data[1]);
      V.v[2] = (double) bswap(data[2]);
#else
      V.v[0] = (double) data[0];
      V.v[1] = (double) data[1];
      V.v[2] = (double) data[2];
#endif

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

    editor.add_cell(c_index++, &index[0]);

    /* Aux data */
    fp.read((char *) &hdr, 2 * sizeof(char));
  }

  editor.init_vertices(vertices.size());
  for (_ordered_set<stl_vertex>::iterator it = vertices.begin();
       it != vertices.end(); ++it)
  {
    editor.add_vertex(it->index, &it->v[0]);
  }

  fp.close();
  editor.close();

}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
