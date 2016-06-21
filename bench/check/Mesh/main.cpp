#include <dolfin.h>

#include <dolfin/common/AdjacentMapping.h>
#include <dolfin/mesh/GlobalFacetMap.h>

using namespace dolfin;

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

    //
    T.begin("Save partitions");
    {
      MeshFunction<uint> p(mesh, mesh.topology().dim());
      p = MPI::processNumber();
      File fp("partitions.pvd");
      fp << p;
    }
    T.end();
    //
    T.begin("Initialize all connectivities");
    {
      mesh.init();
    }
    T.end();
    //
    T.begin("Check");
    {
      check(mesh);
    }
    T.end();
    //
    T.begin("BoundaryMesh exterior");
    {
      BoundaryMesh boundary(mesh, BoundaryMesh::exterior);
      File f("bme.pvd");
      f << boundary;

      message("Initialize all boundary connectivities");
      boundary.init();

      message("Check boundary mesh");
      check(boundary);

      //
      message("Loop on facets");
      uint const tdim = mesh.topology().dim();
      for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
      {
        Facet facet(mesh, boundary.facet_index(*bcell));
        Cell cell(mesh, facet.entities(tdim)[0]);
        uint local_facet = cell.index(facet);

      }
    }
    T.end();
    //
    T.begin("BoundaryMesh interior");
    {
      BoundaryMesh boundary(mesh, BoundaryMesh::interior);
      File f("bmi.pvd");
      f << boundary;

      message("Initialize all boundary connectivities");
      boundary.init();

      message("Check boundary mesh");
      check(boundary);

      //
      message("Loop on facets");
      uint const tdim = mesh.topology().dim();
      for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
      {
        Facet facet(mesh, boundary.facet_index(*bcell));
        Cell cell(mesh, facet.entities(tdim)[0]);
        uint local_facet = cell.index(facet);

      }
    }
    T.end();
  }
  //---------------------------------------------------------------------------
  return 0;
}
