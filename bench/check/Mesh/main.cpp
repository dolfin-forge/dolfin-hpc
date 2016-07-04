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

    /*
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
    //
    T.begin("Compute BoundaryMesh interior/exterior with Hole");
    {
      Hole h;
      BoundaryMesh bint(mesh, h, BoundaryMesh::interior);
      File fi("bmhi.pvd");
      fi << bint;
      BoundaryMesh bext(mesh, h, BoundaryMesh::exterior);
      File fe("bmhe.pvd");
      fe << bext;

      BoundaryMesh bint1(bext, BoundaryMesh::interior);
      File fei("bmhei.pvd");
      fei << bint1;
      BoundaryMesh bext1(bext, BoundaryMesh::exterior);
      File fee("bmhee.pvd");
      fee << bext1;
    }
    T.end();
    //
    T.begin("Compute interior of BoundaryMesh exterior");
    {
      BoundaryMesh boundary0(mesh, BoundaryMesh::exterior);
      BoundaryMesh boundary1(boundary0, BoundaryMesh::interior);
      boundary1.disp();
      File f("bmei.pvd");
      f << boundary1;
    }
    T.end();
    //
    T.begin("Compute Hole subdomain inside of BoundaryMesh exterior");
    {
      BoundaryMesh boundary0(mesh, BoundaryMesh::exterior);
      Hole h;
      BoundaryMesh boundary1(boundary0 ,h, true);
      boundary1.disp();
      File f("bmeshinside.pvd");
      f << boundary1;
    }
    T.end();
    //
    T.begin("Compute Hole subdomain overlap of BoundaryMesh exterior");
    {
      BoundaryMesh boundary0(mesh, BoundaryMesh::exterior);
      Hole h;
      BoundaryMesh boundary1(boundary0, h, false);
      boundary1.disp();
      File f("bmeshoverlap.pvd");
      f << boundary1;
    }
    T.end();
    */
    //
    T.begin("Compute normal of BoundaryMesh exterior");
    {
      uint const rank = MPI::processNumber();
      uint const pe_size = MPI::numProcesses();
      BoundaryMesh boundary(mesh, BoundaryMesh::exterior);
      bool const on_boundary = true;
      SubDomain * subdomain = NULL;

      uint const tdim = mesh.topology().dim();
      uint const gdim = mesh.geometry().dim();

      // Collect facet data and entities inside the subdomain i.e with all
      // vertices contained in the subdomain
      bool * facet_select = new bool[boundary.num_cells()];
      std::fill_n(facet_select, boundary.num_cells(), false);
      real * facet_weight = new real[boundary.num_cells()];
      real * facet_normal = new real[gdim * boundary.num_cells()];
      uint facet_count = 0;

      if (boundary.is_distributed())
      {
        message("Collect shared vertex facets by adjacent ranks");
        _set<uint> const& adjs = boundary.distdata()[0].get_adj_ranks();
        uint * adjranks = new uint[adjs.size()];
        std::copy(adjs.begin(), adjs.end(), adjranks);

        Array<uint> * sendbuf_u = new Array<uint> [pe_size];
        _map<uint, _set<uint> > facet_queue;
        uint * vertex_facets =
            new uint[boundary.topology()(0, tdim - 1).max_connections()];
        for (SharedIterator it(boundary.distdata()[0]); !it.end(); ++it)
        {
          Vertex v(boundary, it.index());
          if ((subdomain != NULL) && (!subdomain->inside(v.x(), on_boundary)))
          {
            continue;
          }
          uint num_facets = 0;
          for (CellIterator bcell(v); !bcell.end(); ++bcell)
          {
            if (!facet_select[bcell->index()] &&
                ((subdomain == NULL) || (subdomain->inside(*bcell, on_boundary))))
            {
              uint const bcell_index = bcell->index();
              facet_select[bcell_index] = true;
              Facet facet(mesh, boundary.facet_index(*bcell));
              Cell cell(mesh, facet.entities(tdim)[0]);
              uint const local_facet = cell.index(facet);
              facet_weight[bcell_index] = cell.facet_area(local_facet);
              cell.normal(local_facet , &facet_normal[gdim * bcell_index]);
              ++facet_count;
              //
              vertex_facets[num_facets] = bcell_index;
              ++num_facets;
            }
          }
          if (num_facets > 0)
          {
            _set<uint> const& vadjs = it.adj();
            for (_set<uint>::const_iterator a = vadjs.begin(); a != vadjs.end();
                 ++a)
            {
              sendbuf_u[*a].push_back(it.global_index());
              sendbuf_u[*a].push_back(num_facets);
              for (uint i = 0; i < num_facets; ++i)
              {
                sendbuf_u[*a].push_back(vertex_facets[i]);
                facet_queue[*a].insert(vertex_facets[i]);
              }
            }
          }
        }
        facet_queue.clear();
        delete [] vertex_facets;
        message("Number of selected %u facets.", facet_count);

        //
        Array<real> * sendbuf_r = new Array<real>[adjs.size()];
        MPI_Status  * status = new MPI_Status[adjs.size()];
        MPI_Request * sendreq_u = new MPI_Request[adjs.size()];
        uint * sendsize_u = new uint[adjs.size()];
        MPI_Request * sendreq_r = new MPI_Request[adjs.size()];
        uint * sendsize_r = new uint[adjs.size()];
        MPI_Request * recvreq_u = new MPI_Request[adjs.size()];
        uint * recvsize_u = new uint[adjs.size()];
        MPI_Request * recvreq_r = new MPI_Request[adjs.size()];
        uint * recvsize_r = new uint[adjs.size()];

        //
        message("Collect facet data by adjacent ranks");
        for (uint i = 0; i < adjs.size(); ++i)
        {
          uint const a = adjranks[i];

          // Pack data by adjacent rank
          // U: foreach vertex [global index, num facets, facets local indices]
          //    + [ all facet local indices ] + number of facets
          // R: foreach facet  [weight, normal components]
          for (_set<uint>::const_iterator it = facet_queue[a].begin();
               it != facet_queue[a].end(); ++it)
          {
            sendbuf_u[a].push_back(*it);
            sendbuf_r[i].push_back(facet_weight[*it]);
            for (uint d = 0; d < gdim; ++d)
            {
              sendbuf_r[i].push_back(facet_normal[gdim*(*it) + d]);
            }
          }
          sendbuf_u[a].push_back(facet_queue[a].size());

          // Exchange buffer sizes
          sendsize_u[i] = sendbuf_u[a].size();
          MPI_Isend(&sendsize_u[i], 1, MPI_UNSIGNED, a, 0, MPI::DOLFIN_COMM,
                    &sendreq_u[i]);
          MPI_Irecv(&recvsize_u[i], 1, MPI_UNSIGNED, a, 0, MPI::DOLFIN_COMM,
                    &recvreq_u[i]);
          sendsize_r[i] = sendbuf_r[i].size();
          MPI_Isend(&sendsize_r[i], 1, MPI_UNSIGNED, a, 1, MPI::DOLFIN_COMM,
                    &sendreq_r[i]);
          MPI_Irecv(&recvsize_r[i], 1, MPI_UNSIGNED, a, 1, MPI::DOLFIN_COMM,
                    &recvreq_r[i]);
        }
        MPI_Waitall(adjs.size(), &sendreq_u[0],&status[0]);
        MPI_Waitall(adjs.size(), &recvreq_u[0],&status[0]);
        MPI_Waitall(adjs.size(), &sendreq_r[0],&status[0]);
        MPI_Waitall(adjs.size(), &recvreq_r[0],&status[0]);

        //
        Array<uint> * recvbuf_u = new Array<uint>[adjs.size()];
        Array<real> * recvbuf_r = new Array<real>[adjs.size()];
        for (uint i = 0; i < adjs.size(); ++i)
        {
          uint const a = adjranks[i];
          MPI_Isend(&sendbuf_u[a][0], sendsize_u[i], MPI_UNSIGNED, a, 0,
                    MPI::DOLFIN_COMM, &sendreq_u[i]);
          recvbuf_u[i].resize(recvsize_u[i]);
          MPI_Irecv(&recvbuf_u[i][0], recvsize_u[i], MPI_UNSIGNED, a, 0,
                    MPI::DOLFIN_COMM, &recvreq_u[i]);
          MPI_Isend(&sendbuf_r[i][0], sendsize_r[i], MPI_DOUBLE, a, 1,
                    MPI::DOLFIN_COMM, &sendreq_r[i]);
          recvbuf_r[i].resize(recvsize_r[i]);
          MPI_Irecv(&recvbuf_r[i][0], recvsize_r[i], MPI_DOUBLE, a, 1,
                    MPI::DOLFIN_COMM, &recvreq_r[i]);
        }

        // Compute vertex normal for inner vertices
        for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
        {
          if (!facet_select[bcell->index()] &&
              ((subdomain == NULL) || (subdomain->inside(*bcell, on_boundary))))
          {
            uint const bcell_index = bcell->index();
            facet_select[bcell_index] = true;
            Facet facet(mesh, boundary.facet_index(*bcell));
            Cell cell(mesh, facet.entities(tdim)[0]);
            uint const local_facet = cell.index(facet);
            facet_weight[bcell_index] = cell.facet_area(local_facet);
            cell.normal(local_facet , &facet_normal[gdim * bcell_index]);
            ++facet_count;
          }
          //
        }
        message("Number of selected %u facets.", facet_count);
        //

        MPI_Waitall(adjs.size(), &sendreq_u[0],&status[0]);
        MPI_Waitall(adjs.size(), &recvreq_u[0],&status[0]);
        MPI_Waitall(adjs.size(), &sendreq_r[0],&status[0]);
        MPI_Waitall(adjs.size(), &recvreq_r[0],&status[0]);

        //
        for (uint i = 0; i < adjs.size(); ++i)
        {
          uint const a = adjranks[i];
          uint const size = recvsize_u[i];
          uint const fdata_size = recvbuf_u[i][size - 1];
          uint const vdata_size = size - 1 - fdata_size;
          dolfin_assert(fdata_size * (gdim + 1) == recvsize_r[i]);
          uint const * facet_data = &recvbuf_u[i][vdata_size];
          // Construct map for facet data
          _map<uint, real *> fmap;
          for (uint k = 0; k < fdata_size; ++k)
          {
            fmap.insert(std::pair<uint, real *>(facet_data[k],
                                                &recvbuf_r[i][k*(gdim + 1)]));
          }
          // Add vertex facets
          uint const * vertex_data = &recvbuf_u[i][0];
          for (uint k = 0; k < vdata_size; k += (2+vertex_data[k + 1]))
          {
            for (uint f = 1; f <= vertex_data[k + 1]; ++f)
            {
              uint const facet_index = vertex_data[k + 1 + f];
              real const * wn = &fmap[facet_index][0];
            }
          }
        }

        delete[] recvbuf_r;
        delete[] recvbuf_u;
        delete[] recvsize_r;
        delete[] recvreq_r;
        delete[] recvsize_u;
        delete[] recvreq_u;
        delete[] sendsize_r;
        delete[] sendreq_r;
        delete[] sendsize_u;
        delete[] sendreq_u;
        delete[] status;

        delete[] sendbuf_r;
        delete[] sendbuf_u;
        delete[] adjranks;

        /*
        Array<uint> * sendbuf_u = new Array<uint> [pe_size];
        Array<real> * sendbuf_r = new Array<real> [pe_size];
        for (GhostIterator it(boundary.distdata()[0]); !it.end(); ++it)
        {
          uint const owner = it.owner();
          Vertex v(boundary, it.index());
          uint facets = 0;
          for (CellIterator c(v); !c.end(); ++c)
          {
            if (facet_select[c->index()])
            {
              ++facets;
              // Collect vertex facets
              sendbuf_r[owner].push_back(facet_weight[it.index()]);
              for (uint d = 0; d < gdim; ++d)
              {
                sendbuf_r[owner].push_back(facet_normal[gdim*it.index() + d]);
              }
            }
          }
          if(facets > 0)
          {
            sendbuf_u[owner].push_back(it.global_index());
            sendbuf_u[owner].push_back(facets);
          }
        }

        //
        _set<uint> const& adjs = boundary.distdata()[0].get_adj_ranks();
        MPI_Status status;
        MPI_Request * request_u = new MPI_Request[adjs.size()];
        MPI_Request * request_r = new MPI_Request[adjs.size()];
        uint recvmax_u = 0;
        uint recvmax_r = 0;
        for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it)
        {
          recvmax_u = std::max(recvmax_u, (uint) sendbuf_u[*it].size());
          recvmax_r = std::max(recvmax_r, (uint) sendbuf_r[*it].size());
        }
        MPI::numGlobalMax(recvmax_u, recvmax_u);
        message("recvmax_u = %u", recvmax_u);
        MPI::numGlobalMax(recvmax_r, recvmax_r);
        message("recvmax_r = %u", recvmax_r);
        uint * recvbuf_u = new uint[recvmax_u];
        real * recvbuf_r = new real[recvmax_r];
        //
        uint i = 0;
        for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it, ++i)
        {
          MPI_Isend(&sendbuf_u[*it][0], sendbuf_u[*it].size(), MPI_UNSIGNED, *it,
                    0, MPI::DOLFIN_COMM, &request_u[i]);
          MPI_Irecv(&recvbuf_u[0], recvmax_u, MPI_UNSIGNED, *it, 0,
                    MPI::DOLFIN_COMM, &request_u[i]);
          MPI_Isend(&sendbuf_r[*it][0], sendbuf_r[*it].size(), MPI_DOUBLE, *it,
                    1, MPI::DOLFIN_COMM, &request_r[i]);
          MPI_Irecv(&recvbuf_r[0], recvmax_r, MPI_UNSIGNED, *it, 1,
                    MPI::DOLFIN_COMM, &request_r[i]);
        }

        //


        // Collect facet data from other ranks
        i = 0;
        for (_set<uint>::const_iterator it = adjs.begin(); it != adjs.end(); ++it, ++i)
        {
          int recvcount;
          MPI_Wait(&request_u[i],&status);
          MPI_Wait(&request_r[i],&status);
          MPI_Get_count(&status, MPI_UNSIGNED, &recvcount);
          real * dataptr = &recvbuf_r[0];
          for (uint k = 0; k < recvcount; k+=2)
          {
            message("recv facet data for %u", recvbuf_u[k]);
            for (uint l = 0; l < recvbuf_u[k + 1]; ++l)
            {
              message("facet weight = %f", *dataptr);
              ++dataptr;
              for (uint d = 0; d < gdim; ++d, ++dataptr)
              {
                message("n_%u = %f", d, *dataptr);
              }
            }
          }
        }

        delete[] recvbuf_r;
        delete[] recvbuf_u;
        delete[] sendbuf_r;
        delete[] sendbuf_u;
        delete[] request_r;
        delete[] request_u;
        */
      }

      delete [] facet_normal;
      delete [] facet_weight;
      delete [] facet_select;

    }
    T.end();
  }
  //---------------------------------------------------------------------------
  return 0;
}
