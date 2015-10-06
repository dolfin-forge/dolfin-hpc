#include <dolfin.h>

using namespace dolfin;

int main(int argc, char** argv)
{
  int ret = 0;

  //---------------------------------------------------------------------------
  Test t(argc, argv);
  Mesh mesh(t.args.mesh_file);
  bool throw_error = true;

  t.begin("Interior boundary");
  {
    BoundaryMesh& boundary = mesh.interior_boundary();

    uint vInvalid = 0;
    for (VertexIterator bvertex(boundary); !bvertex.end(); ++bvertex)
    {
      Vertex v(mesh, boundary.vertex_index(*bvertex));
      if (!v.is_shared())
      {
        ++vInvalid;
      }
    }
    if (vInvalid)
    {
      error("Vertices on the interior boundary are not shared: %d", vInvalid);
    }

    uint const facet_dim = boundary.topology().dim();
    uint const numBoundaryFacets = boundary.numCells();
    uint const numMDSharedFacets = mesh.topology().num_shared(facet_dim);
    if (numBoundaryFacets != numMDSharedFacets)
    {
      error("Boundary facets count not equal to shared facets: %d != %d",
            numBoundaryFacets, numMDSharedFacets);
    }
    uint fInvalid = 0;
    for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
    {
      Facet f(mesh, boundary.facet_index(*bcell));
      if (!f.is_shared())
      {
        ++fInvalid;
      }
    }
    if (fInvalid)
    {
      error("Facets on the interior boundary are not shared: %d", fInvalid);
    }
  }
  t.end();

  t.begin("Exterior boundary");
  {
    BoundaryMesh& boundary = mesh.exterior_boundary();
    uint invalid = 0;
    for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
    {
      Facet f(mesh, boundary.facet_index(*bcell));
      if (f.is_shared())
      {
        ++invalid;
      }
    }
    if (invalid)
    {
      error("Facets on the exterior boundary are shared: %d", invalid);
    }
  }
  t.end();

  return ret;
}

