#include <dolfin.h>

#include <dolfin/common/AdjacentMapping.h>
#include <dolfin/mesh/GlobalFacetMap.h>

using namespace dolfin;

int main(int argc, char** argv)
{
  int ret = 0;

  //---------------------------------------------------------------------------
  Test t(argc, argv);
  {
    Mesh mesh(t.args.mesh_file);
    logm.file();

    if (mesh.is_distributed())
    {

      uint const tdim = mesh.topology().dim();
      for (uint i = 0; i < tdim; ++i)
      {
        message("Check ownership for topological dimension %u", i);
        mesh.init(i);

        DistributedData& distdata = mesh.distdata()[i];

        //
        message("Check owner");
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
        message("recvmax %u", recvmax);
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

          uint not_owned = 0;
          for (int k = 0; k < recvcount; ++k)
          {
            //
            if(!distdata.has_global(recvbuf[k]))
            {
              ++not_owned;
            }

            //
            uint const local_index = distdata.get_local(recvbuf[k]);
            if(distdata.get_shared_adj(local_index).count(src) == 0)
            {
              error("Missing adjacency");
            }
          }

          if(not_owned > 0)
          {
            warning("Found %u  entities not owned by %u", not_owned, rank);
          }
        }
        delete[] recvbuf;
        delete[] sendbuf;

        //
        message("SharedMapping");
        SharedMapping sm(distdata);
        sm.disp();

        //
        _set<uint> const& adj_ranks = distdata.get_adj_ranks();

      }

      //
      message("GlobalFacetMap");
      GlobalFacetMap gfm(mesh);
      gfm.disp();

      //
      //message("NodeNormal");
      //NodeNormal nn(mesh);
    }
  }
  //---------------------------------------------------------------------------
  return 0;
}
