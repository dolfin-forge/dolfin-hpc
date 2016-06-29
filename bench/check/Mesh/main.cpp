#include <dolfin.h>

#include <dolfin/common/AdjacentMapping.h>
#include <dolfin/mesh/GlobalFacetMap.h>
#include <dolfin/mesh/SubDomain.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
class Hole : public SubDomain
{

public:

  bool inside(real const * x, bool const on_boundary) const
  {

    return (on_boundary
        && ((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.0) * (x[1] - 0.0)) > 0.1);
  }

};
//-----------------------------------------------------------------------------
void check(Mesh& mesh)
{
  if (mesh.is_distributed())
  {

    uint const rank = MPI::processNumber();
    uint const pe_size = MPI::numProcesses();
    uint const tdim = mesh.topology().dim();
    for (uint i = 0; i < tdim; ++i)
    {
      message("Initialize topological dimension %u", i);
      mesh.init(i);

      DistributedData& distdata = mesh.distdata()[i];

      //
      message("Check ownership for topological dimension %u", i);
      {
        uint const rank = MPI::processNumber();
        uint const pe_size = MPI::numProcesses();
        Array<uint> * sendbuf = new Array<uint> [pe_size];
        for (GhostIterator it(distdata); !it.end(); ++it)
        {
          sendbuf[it.owner()].push_back(it.global_index());
        }

        //
        MPI_Status status;
        int src;
        int dst;
        uint recvmax = sendbuf[0].size();
        for (uint j = 0; j < pe_size; ++j)
        {
          recvmax = std::max((uint) sendbuf[j].size(), recvmax);
        }
        MPI::numGlobalMax(recvmax, recvmax);
        uint * recvbuf = new uint[recvmax];
        int recvcount;
        for (uint j = 1; j < pe_size; ++j)
        {
          src = (rank - j + pe_size) % pe_size;
          dst = (rank + j) % pe_size;

          MPI_Sendrecv(&sendbuf[dst][0], sendbuf[dst].size(), MPI_UNSIGNED, dst,
                       0, &recvbuf[0], recvmax, MPI_UNSIGNED, src, 0,
                       MPI::DOLFIN_COMM, &status);
          MPI_Get_count(&status, MPI_UNSIGNED, &recvcount);

          uint not_known = 0;
          uint not_owned = 0;
          for (int k = 0; k < recvcount; ++k)
          {
            //
            if (!distdata.has_global(recvbuf[k]))
            {
              ++not_known;
              continue;
            }
            //
            uint const local_index = distdata.get_local(recvbuf[k]);
            if (!distdata.is_owned(local_index))
            {
              ++not_owned;
            }
            //
            if (distdata.get_shared_adj(local_index).count(src) == 0)
            {
              error("Missing adjacency");
            }
          }

          if (not_owned > 0)
          {
            error("Found %u  entities not owned by %u", not_owned, rank);
          }
        }
        delete[] recvbuf;
        delete[] sendbuf;
      }
      //
      message("Test adjacent ranks for topological dimension %u", i);
      {

        _set<uint> const& adjs = distdata.get_adj_ranks();
        message("On rank %u", rank);
        for (_set<uint>::const_iterator adj = adjs.begin(); adj != adjs.end();
            ++adj)
        {
          if (*adj == rank)
          {
            error("Self is adjacent of self");
          }
          message("%u", *adj);
        }
        Array<uint> sendbuf;
        sendbuf.assign(adjs.begin(), adjs.end());
        MPI_Status status;
        uint recvmax = sendbuf.size();
        MPI::numGlobalMax(recvmax,recvmax);
        uint * recvbuf = new uint[recvmax];

        for (_set<uint>::const_iterator adj = adjs.begin(); adj != adjs.end();
            ++adj)
        {
          MPI_Send(&sendbuf[0], sendbuf.size(), MPI_UNSIGNED, (*adj), 0,
              MPI::DOLFIN_COMM);
        }
        for (_set<uint>::const_iterator adj = adjs.begin(); adj != adjs.end();
            ++adj)
        {
          MPI_Recv(&recvbuf[0], recvmax, MPI_UNSIGNED, (*adj), 0,
              MPI::DOLFIN_COMM, &status);
          int recvcount;
          MPI_Get_count(&status, MPI_UNSIGNED, &recvcount);
          bool found = false;
          message("Receiving from %u", *adj);
          for (uint k = 0; k < recvcount; ++k)
          {
            message("%u\n", recvbuf[k]);
            if (recvbuf[k] == rank)
            {
              found = true;
              break;
            }
          }
          if (!found)
          {
            error("Adjacency: current rank %u not on src %u", rank, *adj);
          }
        }

        delete [] recvbuf;
      }
      //
      message("Test SharedMapping for topological dimension %u", i);
      {
        SharedMapping sm(distdata);
      }
      //
    }
  }
}
//-----------------------------------------------------------------------------

int main(int argc, char** argv)
{
  int ret = 0;

  //---------------------------------------------------------------------------
  Test T(argc, argv);
  {
    Mesh mesh(T.args.mesh_file);
    logm.verbose(1);
    logm.file();

//    //
//    T.begin("Save partitions");
//    {
//      MeshFunction<uint> p(mesh, mesh.topology().dim());
//      p = MPI::processNumber();
//      File fp("partitions.pvd");
//      fp << p;
//    }
//    T.end();
//    //
//    T.begin("Initialize all connectivities");
//    {
//      mesh.init();
//    }
//    T.end();
//    //
//    T.begin("Check");
//    {
//      check(mesh);
//    }
//    T.end();
//    //
//    T.begin("BoundaryMesh exterior");
//    {
//      BoundaryMesh boundary(mesh, BoundaryMesh::exterior);
//      File f("bme.pvd");
//      f << boundary;
//
//      message("Initialize all boundary connectivities");
//      boundary.init();
//
//      message("Check boundary mesh");
//      check(boundary);
//
//      //
//      message("Loop on facets");
//      uint const tdim = mesh.topology().dim();
//      for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
//      {
//        Facet facet(mesh, boundary.facet_index(*bcell));
//        Cell cell(mesh, facet.entities(tdim)[0]);
//        uint local_facet = cell.index(facet);
//
//      }
//    }
//    T.end();
//    //
//    T.begin("BoundaryMesh interior");
//    {
//      BoundaryMesh boundary(mesh, BoundaryMesh::interior);
//      File f("bmi.pvd");
//      f << boundary;
//
//      message("Initialize all boundary connectivities");
//      boundary.init();
//
//      message("Check boundary mesh");
//      check(boundary);
//
//      //
//      message("Loop on facets");
//      uint const tdim = mesh.topology().dim();
//      for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
//      {
//        Facet facet(mesh, boundary.facet_index(*bcell));
//        Cell cell(mesh, facet.entities(tdim)[0]);
//        uint local_facet = cell.index(facet);
//
//      }
//    }
//    T.end();
//    //
//    T.begin("Compute BoundaryMesh interior/exterior with Hole");
//    {
//      Hole h;
//      BoundaryMesh bint(mesh, h, BoundaryMesh::interior);
//      File fi("bmhi.pvd");
//      fi << bint;
//      BoundaryMesh bext(mesh, h, BoundaryMesh::exterior);
//      File fe("bmhe.pvd");
//      fe << bext;
//
//      BoundaryMesh bint1(bext, BoundaryMesh::interior);
//      File fei("bmhei.pvd");
//      fei << bint1;
//      BoundaryMesh bext1(bext, BoundaryMesh::exterior);
//      File fee("bmhee.pvd");
//      fee << bext1;
//    }
//    T.end();
//    //
//    T.begin("Compute interior of BoundaryMesh exterior");
//    {
//      BoundaryMesh boundary0(mesh, BoundaryMesh::exterior);
//      BoundaryMesh boundary1(boundary0, BoundaryMesh::interior);
//      boundary1.disp();
//      File f("bmei.pvd");
//      f << boundary1;
//    }
//    T.end();
    //
    T.begin("Compute normal of BoundaryMesh exterior");
    {
      uint const rank = MPI::processNumber();
      uint const pe_size = MPI::numProcesses();
      BoundaryMesh boundary(mesh, BoundaryMesh::exterior);

      uint const tdim = mesh.topology().dim();
      uint const gdim = mesh.geometry().dim();

      //
      real * facet_weights = new real[boundary.num_cells()];
      real * facet_normals = new real[gdim * boundary.num_cells()];
      for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
      {
        Facet facet(mesh, boundary.facet_index(*bcell));
        Cell cell(mesh, facet.entities(tdim)[0]);
        uint const local_facet = cell.index(facet);
        facet_weights[bcell->index()] = cell.facet_area(local_facet);
        cell.normal(bcell->index())
        std::copy(  facet_weights[gdim*bcell->index()]);
        //
        if(cell.entities(tdim - 1)[local_facet] != boundary.facet_index(*bcell))
        {
          error("%u != %u", cell.entities(tdim - 1)[local_facet], boundary.facet_index(*bcell));
        }
        //
        if(!abscmp(bcell->volume(), cell.facet_area(local_facet)))
        {
          message("%f == %f", bcell->volume(), cell.facet_area(local_facet));
          error("Inconsistent facet computation");
        }
        //
        for (VertexIterator v(*bcell); !v.end(); ++v)
        {
          message("%u:", boundary.vertex_index(v->index()));
          v->point().disp();
        }
        for (VertexIterator v(facet); !v.end(); ++v)
        {
          message("%u", v->index());
          v->point().disp();
        }
        skip();

      }
      delete[] facet_normals;
      delete[] facet_weights;

      if (mesh.is_distributed())
      {
        message("Collect facets per vertices");
        Array<uint> * sendbuf_u = new Array<uint> [pe_size];
        for (GhostIterator it(boundary.distdata()[0]); !it.end(); ++it)
        {
          sendbuf_u[it.owner()].push_back(it.global_index());
          //

        }
        delete[] sendbuf_u;

        /*
         for (uint d = 0; d <= tdim - 1; ++d)
         {
         // Initialize topological dimension
         boundary.init(d);
         Array<uint> * ghost_entities = new Array<uint> [pe_size];
         for (GhostIterator it(boundary.distdata()[d]); !it.end(); ++it)
         {

         }

         for (MeshEntityIterator it(boundary, d); !it.end(); ++it)
         {
         if (it->is_owned())
         {

         }
         else
         {
         //            ghost_entities[it->owner()].push_back()
         }
         }
         delete[] ghost_entities;
         }
         */

      }

    }
    T.end();
  }
  //---------------------------------------------------------------------------
  return 0;
}
