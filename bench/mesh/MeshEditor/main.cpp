#include <dolfin.h>

using namespace dolfin;

int main(int argc, char** argv)
{
  dolfin_init(argc, argv);
  //---------------------------------------------------------------------------
  Mesh mesh;
  Array<CellType::Type> types;
  types.push_back(CellType::point);
  types.push_back(CellType::interval);
  types.push_back(CellType::triangle);
  types.push_back(CellType::tetrahedron);
  uint const N = 16;

  for (Array<CellType::Type>::const_iterator it = types.begin();
      it != types.end(); ++it)
  {
    CellType * celltype = CellType::create(*it);
    uint const tdim = celltype->dim();
    for (uint gdim = tdim; gdim <= EuclideanSpace::MAX_DIMENSION; ++gdim)
    {
      message("Create mesh editor:\n\t'%s' in R^%d",
              celltype->description().c_str(), gdim);
      begin("");

      //---
      tic();
      MeshEditor editor(mesh, celltype->cellType(), gdim);

      uint const nvertices = (tdim == 0 ? 1 : std::pow(N + 1, tdim));
      message("Initialize vertices : %d", nvertices);
      editor.init_vertices(nvertices);

      uint const ncells = (tdim == 0 ? 1 : fact(tdim) * std::pow(N, tdim));
      message("Initialize cells    : %d", ncells);
      editor.init_cells(ncells);

      // Vertex coordinates
      real * coords = new real[nvertices*gdim];
      uint x = 0;
      message("Adding coordinates in R^%d.", gdim);
      for (uint v = 0; v < nvertices; ++v)
      {
        for(uint d = 0; d < gdim; ++d)
        {
          coords[x++] = (real) d;
        }
      }
      uint vertexid = 0;
      message("Adding vertices.");
      for (uint v = 0; v < nvertices; ++v)
      {
        editor.add_vertex(vertexid++, &coords[v * gdim]);
      }
      delete[] coords;

      // Connectivities
      uint cellid = 0;
      uint * connectivities = new uint[celltype->num_entities(0)];
      message("Adding connectivities: %d per cell.", celltype->num_entities(0));
      for (uint c = 0; c < ncells; ++c)
      {
        editor.add_cell(cellid++, &connectivities[0]);
      }
      delete[] connectivities;

      editor.close();
      tocd();
      //---

      message("Done");
      end();
      skip();
    }

    delete celltype;
    celltype = NULL;
  }

  dolfin_finalize();
  return 0;
}

