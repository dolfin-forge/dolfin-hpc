
#include <dolfin.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
inline real length_old(Edge& e)
{
  uint* vertices = e.entities(0);
  dolfin_assert(vertices);

  const Vertex v0(e.mesh(), vertices[0]);
  const Vertex v1(e.mesh(), vertices[1]);

  const Point p0 = v0.point();
  const Point p1 = v1.point();

  real length(
      sqrt(
          (p1.x() - p0.x()) * (p1.x() - p0.x())
              + (p1.y() - p0.y()) * (p1.y() - p0.y())
              + (p1.z() - p0.z()) * (p1.z() - p0.z())));

  return length;
}
//-----------------------------------------------------------------------------
inline real length_loop(Edge const& e)
{
  uint const * vertices = e.entities(0);
  MeshGeometry const& geom = e.mesh().geometry();
  uint const d = geom.dim();
  real const * p0 = geom.x(vertices[0]);
  real const * p1 = geom.x(vertices[1]);
  real l = 0;
  for (uint n = 0; n < d; ++n)
  {
    l += (p1[n] - p0[n]) * (p1[n] - p0[n]);
  }
  return std::sqrt(l);
}
//-----------------------------------------------------------------------------
inline real length_memcpy(Edge const& e)
{
  MeshGeometry const& geom = e.mesh().geometry();
  uint const d = geom.dim();
  uint const * vertices = e.entities(0);
  real p0[3] = { 0.0 };
  real p1[3] = { 0.0 };
  memcpy(&p0, geom.x(vertices[0]), d * sizeof(real));
  memcpy(&p1, geom.x(vertices[1]), d * sizeof(real));
  return std::sqrt(
      (p1[0] - p0[0]) * (p1[0] - p0[0]) + (p1[1] - p0[1]) * (p1[1] - p0[1])
          + (p1[2] - p0[2]) * (p1[2] - p0[2]));
}
//-----------------------------------------------------------------------------
Point midpoint_old(Edge& e)
{
  uint* vertices = e.entities(0);
  dolfin_assert(vertices);

  const Vertex v0(e.mesh(), vertices[0]);
  const Vertex v1(e.mesh(), vertices[1]);

  const Point p0 = v0.point();
  const Point p1 = v1.point();

  Point p(0.5 * (p0.x() + p1.x()), 0.5 * (p0.y() + p1.y()),
          0.5 * (p0.z() + p1.z()));

  return p;
}

int main(int argc, char** argv)
{
  Mesh mesh("square4M.bin");

  message("Initialize vertex-edge connectivities");
  mesh.init(0, 1);

  message("Number of cells in the mesh: %d.", mesh.num_global_cells());
  message("Number of edges in the mesh: %d.", mesh.global_size(1));

  uint const N = 1e2;
  message("Repeating %d times.", N);

  message("Edge: create first entity");
  Edge e0(mesh, 0);

  message("EdgeIterator: loop on mesh entities");
  tic();
  for (uint i = 0; i < N; ++i)
  {
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {

    }
  }
  tocd();

  //---------------------------------------------------------------------------

  real * Alength = new real[mesh.size(1)];

  message("EdgeIterator: compute length (class)");
  tic();
  for (uint i = 0; i < N; ++i)
  {
    uint k = 0;
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      Alength[k++] = e->length();
    }
  }
  tocd();

  message("EdgeIterator: compute length (old function)");
  tic();
  for (uint i = 0; i < N; ++i)
  {
    uint k = 0;
    uint const d = mesh.geometry().dim();
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      real l = length_old(*e);
      dolfin_assert(abscmp(Alength[k++], l));
    }
  }
  tocd();

  message("EdgeIterator: compute length (manual loop)");
  tic();
  for (uint i = 0; i < N; ++i)
  {
    uint k = 0;
    uint const d = mesh.geometry().dim();
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      uint * vertices = e->entities(0);
      real * p0 = mesh.geometry().x(vertices[0]);
      real * p1 = mesh.geometry().x(vertices[1]);
      real l = 0;
      for (uint n = 0; n < d; ++n)
      {
        l += (p1[n] - p0[n]) * (p1[n] - p0[n]);
      }
      l = std::sqrt(l);
      dolfin_assert(abscmp(Alength[k++], l));
    }
  }
  tocd();

  message("EdgeIterator: compute length (function loop)");
  tic();
  for (uint i = 0; i < N; ++i)
  {
    uint k = 0;
    uint const d = mesh.geometry().dim();
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      real l = length_loop(*e);
      dolfin_assert(abscmp(Alength[k++], l));
    }
  }
  tocd();

  message("EdgeIterator: compute length (manual memcpy)");
  tic();
  for (uint i = 0; i < N; ++i)
  {
    uint k = 0;
    uint const d = mesh.geometry().dim();
    real p0[3] = { 0.0 };
    real p1[3] = { 0.0 };
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      uint * vertices = e->entities(0);
      memcpy(&p0, mesh.geometry().x(vertices[0]), d * sizeof(real));
      memcpy(&p1, mesh.geometry().x(vertices[1]), d * sizeof(real));
      real l = std::sqrt(
          (p1[0] - p0[0]) * (p1[0] - p0[0]) + (p1[1] - p0[1]) * (p1[1] - p0[1])
              + (p1[2] - p0[2]) * (p1[2] - p0[2]));
      dolfin_assert(abscmp(Alength[k++], l));
    }
  }
  tocd();

  message("EdgeIterator: compute length (function memcpy)");
  tic();
  for (uint i = 0; i < N; ++i)
  {
    uint k = 0;
    uint const d = mesh.geometry().dim();
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      real l = length_memcpy(*e);
      dolfin_assert(abscmp(Alength[k++], l));
    }
  }
  tocd();

  delete[] Alength;

  //---------------------------------------------------------------------------

  Point * Amidpnt = new Point[mesh.size(1)];

  message("EdgeIterator: compute midpoint (class)");
  tic();
  for (uint i = 0; i < N; ++i)
  {
    uint k = 0;
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      Amidpnt[k++] = e->midpoint();
    }
  }
  tocd();

  message("EdgeIterator: compute midpoint (old function)");
  tic();
  for (uint i = 0; i < N; ++i)
  {
    uint k = 0;
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      Point p = midpoint_old(*e);
      dolfin_assert(p.distance(Amidpnt[k++]) < DOLFIN_EPS);
    }
  }
  tocd();

  message("EdgeIterator: compute midpoint (manual loop)");
  tic();
  for (uint i = 0; i < N; ++i)
  {
    uint k = 0;
    MeshGeometry const& geom = mesh.geometry();
    uint const d = mesh.geometry().dim();
    for (EdgeIterator e(mesh); !e.end(); ++e)
    {
      uint* vertices = e->entities(0);
      dolfin_assert(vertices);
      real const * p0 = geom.x(vertices[0]);
      real const * p1 = geom.x(vertices[1]);
      Point p;
      for (uint n = 0; n < d; ++n)
      {
        p[n] = 0.5 * (p0[n] + p1[n]);
      } dolfin_assert(p.distance(Amidpnt[k++]) < DOLFIN_EPS);
    }
  }
  tocd();

  delete[] Amidpnt;

  return 0;
}

