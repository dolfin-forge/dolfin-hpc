#include <dolfin.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  int ret = 0;

  //---------------------------------------------------------------------------
  Test t(argc, argv);
  Mesh mesh(t.args.mesh_file);
  bool throw_error = true;

  t.begin("MPIMeshCommunicator::distribute(Mesh<uint,Cell>)");
  {
    MeshData D(mesh);
    uint const pe_size = PE::size();
    uint const pe_rank = PE::rank();

    // Cell-based distribution
    MeshValues<uint, Cell> dist(mesh);
    for (Cell::iterator c(mesh); !c.end(); ++c)
    {
      dist(*c) = std::rand() % pe_size;
    }
    D.add(dist);

    // Vertex coordinates
    uint const gdim = mesh.geometry().dim();
    Array< MeshValues<real, Vertex> > x(gdim, MeshValues<real, Vertex>(mesh));
    for (Vertex::iterator v(mesh); !v.end(); ++v)
    {
      for (uint d = 0; d < gdim; ++d) x[d](*v) = v->x()[d];
    }
    for (uint d = 0; d < gdim; ++d) D.add(x[d]);

    // Redistribute according to cell distribution
    mesh.distribute(dist, D);

    // Check transferred MeshValues
    uint invalid;
    {
      invalid= 0;
      for (Cell::iterator c(mesh); !c.end(); ++c)
      {
        if (dist(*c) != pe_rank) invalid++;
      }
      if (invalid) error("Invalid cell transfer: %u", invalid);
    }

    {
      invalid= 0;
      for (Vertex::iterator v(mesh); !v.end(); ++v)
      {
        // If MeshValues are now synchronized then ghosts hold zero values
        if (v->is_owned())
        {
          for (uint d = 0; d < gdim; ++d)
          {
            if (x[d](*v) != v->x()[d])
            {
              invalid++;
            }
          }
        }
      }
      if (invalid) error("Invalid vertex coordinates transfer: %u", invalid);
    }
  }
  t.end();

  return ret;
}

