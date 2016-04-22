#include <dolfin.h>

#include <dolfin/fem/ScratchSpace.h>

using namespace dolfin;



void test_boundary_facets(Mesh& mesh)
{
  BoundaryMesh& exterior = mesh.exterior_boundary();
  uint const rank = dolfin::MPI::processNumber();
  uint const pe_size = dolfin::MPI::numProcesses();
  uint * sendbuf = new uint[pe_size];
  std::fill_n(sendbuf, pe_size, 0);
  _set<uint> const& adjs = exterior.distdata()[0].get_adj_ranks();
  for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it)
  {
    sendbuf[*it] = exterior.num_cells();
  }

  MPI_Status status;
  uint src;
  uint dst;

  uint * recvbuf = new uint[pe_size];
  for (uint j = 1; j < pe_size; ++j)
  {
    src = (rank - j + pe_size) % pe_size;
    dst = (rank + j) % pe_size;

    MPI_Sendrecv(&sendbuf[dst], 1, MPI_UNSIGNED, dst, 1, &recvbuf[src], 1,
    MPI_UNSIGNED,
                 src, 1, dolfin::MPI::DOLFIN_COMM, &status);

  }

  for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it)
  {
    if (sendbuf[*it] == 0)
    {
      error("Adjacent rank %u has no boundary cell", *it);
    }
    else
    {
      message("Adjacent rank %u has %u boundary cell", *it, sendbuf[*it]);
    }

  }

  delete[] recvbuf;
  delete[] sendbuf;

}
//
uint const DEGMIN = 2;
uint const DEGMAX = 2;
void test(std::string name, NodeNormal& vn)
{
  Mesh& mesh = vn.mesh();
  test_boundary_facets(mesh);
  for (uint i = DEGMIN; i <= DEGMAX; ++i)
  {
    message("Computing for CG%u", i);
    uint const tdim = mesh.topology().dim();
    uint const gdim = mesh.geometry().dim();
    ufl::VectorElement cgd(ufl::Family::CG, vn.mesh().type(), i, gdim);
    FiniteElementSpace Vh(vn.mesh(), cgd);

    BoundaryMesh& boundary = mesh.exterior_boundary();
    for (uint d = 1; d < tdim - 1; ++d)
    {
      if (Vh.dofmap().needs_mesh_entities(d))
      {
        boundary.init(d);
        message("boundary entities of dimension %u", d);
        for (SharedIterator it(boundary.distdata()[d]); !it.end(); ++it)
        {
          message("shared index %8u", it.index());
        }
      }
    }

    ScratchSpace scratch(Vh);
    DofMap const& dofmap = Vh.dofmap();
    uint * entity_dofs = new uint[dofmap.local_dimension()];
    for (CellIterator c(boundary); !c.end(); ++c)
    {
      Facet facet(mesh, boundary.facet_index(c->index()));
      Cell cell(mesh, facet.entities(tdim)[0]);
      scratch.cell.update(cell);

      message("cell %u", cell.index());
      for (uint d = 0; d < tdim; ++d)
      {
        message("%u:", d);
        begin("");
        for (MeshEntityIterator it(cell, d); !it.end(); ++it)
        {
          message("%u: %8u", it.pos(), it->index());
        }
        end();
      }

      message("per entity");
      uint local_facet = cell.index(facet);
      dofmap.tabulate_dofs(scratch.dofs, scratch.cell);
      dofmap.tabulate_facet_dofs(scratch.facet_dofs, local_facet);
      for (uint d = 0; d < tdim - 1; ++d)
      {
        if (dofmap.needs_mesh_entities(d))
        {
          message("%u:", d);
          begin("");
          for (MeshEntityIterator it(facet, d); !it.end(); ++it)
          {
            uint const local_index = cell.index(*it);
            message("%u: %8u", local_index, it->index());
            dofmap.tabulate_entity_dofs(entity_dofs, d, local_index);
            if (it->is_ghost())
            {
              for (uint n = 0; n < dofmap.num_entity_dofs(d); ++n)
              {
                if (!dofmap.is_ghost(scratch.dofs[entity_dofs[n]]))
                {
                  error("dofs %u is not ghost", scratch.dofs[entity_dofs[n]]);
                }
                if (!dofmap.is_shared(scratch.dofs[entity_dofs[n]]))
                {
                  error("ghost dofs %u is not shared",
                        scratch.dofs[entity_dofs[n]]);
                }
              }
            }
          }
          end();
        }
      }
    }
    delete[] entity_dofs;

    vn.init(Vh);
    vn.compute();
    std::stringstream ss;
    ss << "CG" << i << name << ".pvd";
    vn.write(ss.str());
  }
}
//
int main(int argc, char** argv)
{
  dolfin_init(argc, argv);
  logm.verbose(1);
  logm.file();
//---------------------------------------------------------------------------
  {
    {
      Mesh mesh("../../data/meshes/squareN100R.xml.gz");
      NodeNormal vn(mesh, NodeNormal::none);
      test("square_none", vn);
    }

    {
      Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
      NodeNormal vn(mesh, NodeNormal::none);
      test("cube_none", vn);
    }

    {
      Mesh mesh("../../data/meshes/squareN100R.xml.gz");
      NodeNormal vn(mesh, NodeNormal::facet);
      test("square_facet", vn);
    }

    {
      Mesh mesh("../../data/meshes/cubeN32R.xml.gz");
      NodeNormal vn(mesh, NodeNormal::facet);
      test("cube_facet", vn);
    }

  }
//---------------------------------------------------------------------------
  dolfin_finalize();
  return 0;
}

