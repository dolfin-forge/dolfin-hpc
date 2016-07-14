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
      bool * facet_computed = new bool[boundary.num_cells()];
      std::fill_n(facet_computed, boundary.num_cells(), false);
      uint const facet_data_size = gdim + 1;
      real * facet_data = new real[facet_data_size * boundary.num_cells()];
      uint facet_count = 0;

      if (boundary.is_distributed())
      {
        _set<uint> const& adjs = boundary.distdata()[0].get_adj_ranks();
        uint * adjranks = new uint[adjs.size()];
        std::copy(adjs.begin(), adjs.end(), adjranks);

        // Collect shared vertex facets by adjacent ranks and order sending
        message("Collect shared vertex facets by adjacent ranks and order sending");
        tic();
        Array<uint> * sendbuf_u = new Array<uint> [pe_size];
        _map<uint, uint> * facet_sendmap = new _map<uint, uint>[pe_size];
        uint * vfacets = new uint[boundary.topology()(0, tdim-1).max_connections()];
        for (SharedIterator it(boundary.distdata()[0]); !it.end(); ++it)
        {
          Vertex v(boundary, it.index());
          if ((subdomain != NULL) && !subdomain->inside(v.x(), on_boundary))
          {
            continue;
          }
          uint num_facets = 0;
          for (CellIterator bcell(v); !bcell.end(); ++bcell)
          {
            if ((subdomain == NULL) || subdomain->inside(*bcell, on_boundary))
            {
              uint const id = bcell->index();
              if(!facet_computed[bcell->index()])
              {
                facet_computed[id] = true;
                Facet facet(mesh, boundary.facet_index(*bcell));
                Cell cell(mesh, facet.entities(tdim)[0]);
                uint const local_facet = cell.index(facet);
                facet_data[facet_data_size * id] = cell.facet_area(local_facet);
                cell.normal(local_facet , &facet_data[facet_data_size * id + 1]);
                ++facet_count;
              }
              //
              vfacets[num_facets++] = id;
            }
          }
          if (num_facets > 0)
          {
            //message("v%8u: num facets = %8u", it.global_index(), num_facets);
            _set<uint> const& vadjs = it.adj();
            for (_set<uint>::const_iterator a = vadjs.begin(); a != vadjs.end();
                 ++a)
            {
              sendbuf_u[*a].push_back(it.global_index());
              sendbuf_u[*a].push_back(num_facets);
              //message("\tadj%8u:", *a);
              for (uint i = 0; i < num_facets; ++i)
              {
                uint const id = vfacets[i];
                _map<uint, uint>::const_iterator it = facet_sendmap[*a].find(id);
                if (it == facet_sendmap[*a].end())
                {
                  uint const count = facet_sendmap[*a].size();
                  sendbuf_u[*a].push_back(count);
                  facet_sendmap[*a][id] = count;
                  //message("\t\tadd facet %8u as %8u (new)", id, count);
                }
                else
                {
                  sendbuf_u[*a].push_back(it->second);
                  //message("\t\tadd facet %8u as %8u (use)", id, it->second);
                }
              }
            }
          }
        }
        delete [] vfacets;
        tocd();

        // Create MPI variables
        message("Create MPI variables");
        MPI_Status  * status = new MPI_Status[adjs.size()];
        MPI_Request * sendreq_u = new MPI_Request[adjs.size()];
        uint * sendsize_u = new uint[adjs.size()];
        MPI_Request * recvreq_u = new MPI_Request[adjs.size()];
        uint * recvsize_u = new uint[adjs.size()];
        MPI_Request * sendreq_r = new MPI_Request[adjs.size()];
        uint * sendsize_r = new uint[adjs.size()];
        MPI_Request * recvreq_r = new MPI_Request[adjs.size()];
        uint * recvsize_r = new uint[adjs.size()];

        // Exchange buffer sizes
        message("Exchange buffer sizes");
        for (uint i = 0; i < adjs.size(); ++i)
        {
          uint const a = adjranks[i];

          // Exchange buffer sizes
          sendsize_u[i] = sendbuf_u[a].size();
          MPI_Isend(&sendsize_u[i], 1, MPI_UNSIGNED, a, 0, MPI::DOLFIN_COMM,
                    &sendreq_u[i]);
          MPI_Irecv(&recvsize_u[i], 1, MPI_UNSIGNED, a, 0, MPI::DOLFIN_COMM,
                    &recvreq_u[i]);

          sendsize_r[i] = facet_sendmap[a].size() * facet_data_size;
          MPI_Isend(&sendsize_r[i], 1, MPI_UNSIGNED, a, 1, MPI::DOLFIN_COMM,
                    &sendreq_r[i]);
          MPI_Irecv(&recvsize_r[i], 1, MPI_UNSIGNED, a, 1, MPI::DOLFIN_COMM,
                    &recvreq_r[i]);
        }
        MPI_Waitall(adjs.size(), &sendreq_u[0],&status[0]);
        MPI_Waitall(adjs.size(), &recvreq_u[0],&status[0]);

        // Start exchange of vertex - facets connectivities
        // U: foreach vertex [global index, num facets, facets send indices]
        message("Start exchange of vertex - facets connectivities");
        uint * offsets_u = new uint[adjs.size() + 1];
        offsets_u[0] = 0;
        for (uint i = 0; i < adjs.size(); ++i)
        {
          offsets_u[i + 1] = offsets_u[i] + recvsize_u[i];
        }
        uint * recvbuf_u = new uint[offsets_u[adjs.size()]];
        for (uint i = 0; i < adjs.size(); ++i)
        {
          uint const a = adjranks[i];
          MPI_Isend(&sendbuf_u[a][0], sendsize_u[i], MPI_UNSIGNED, a, 0,
                    MPI::DOLFIN_COMM, &sendreq_u[i]);
          MPI_Irecv(&recvbuf_u[offsets_u[i]], recvsize_u[i], MPI_UNSIGNED, a, 0,
                    MPI::DOLFIN_COMM, &recvreq_u[i]);
        }

        // Collect facet data by adjacent ranks in order of addition
        message("Collect facet data by adjacent ranks in order of addition");
        tic();
        Array<real> * sendbuf_r = new Array<real>[adjs.size()];
        for (uint i = 0; i < adjs.size(); ++i)
        {
          uint const a = adjranks[i];

          // Pack data by adjacent rank
          // R: foreach facet  [weight, normal components]
          sendbuf_r[i].resize(facet_sendmap[a].size() * facet_data_size);
          dolfin_assert(sendsize_r[i] == sendbuf_r[i].size());
          for (_map<uint, uint>::const_iterator it = facet_sendmap[a].begin();
               it != facet_sendmap[a].end(); ++it)
          {
            dolfin_assert(it->second < sendsize_u[i]);
            std::copy(&facet_data[facet_data_size*it->first],
                      &facet_data[facet_data_size*it->first] + facet_data_size,
                      &sendbuf_r[i][it->second*facet_data_size]);
          }
        }
        delete [] facet_sendmap;
        tocd();

        // Wait for reception of facet data buffer sizes
        message("Wait for reception of facet data buffer sizes");
        MPI_Waitall(adjs.size(), &sendreq_r[0],&status[0]);
        MPI_Waitall(adjs.size(), &recvreq_r[0],&status[0]);

        // Exchange facet data
        message("Exchange facet data");
        uint * offsets_r = new uint[adjs.size() + 1];
        offsets_r[0] = 0;
        for (uint i = 0; i < adjs.size(); ++i)
        {
          offsets_r[i + 1] = offsets_r[i] + recvsize_r[i];
        }
        real * recvbuf_r = new real[offsets_r[adjs.size()]];
        for (uint i = 0; i < adjs.size(); ++i)
        {
          uint const a = adjranks[i];
          MPI_Isend(&sendbuf_r[i][0], sendsize_r[i], MPI_DOUBLE, a, 1,
                    MPI::DOLFIN_COMM, &sendreq_r[i]);
          MPI_Irecv(&recvbuf_r[offsets_r[i]], recvsize_r[i], MPI_DOUBLE, a, 1,
                    MPI::DOLFIN_COMM, &recvreq_r[i]);
        }

        // Compute vertex normal for inner vertices
        message("Compute vertex normal for inner vertices");
        for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
        {
          // Do not recompute weight and normal
          if (!facet_computed[bcell->index()] &&
              ((subdomain == NULL) || subdomain->inside(*bcell, on_boundary)))
          {
            uint const id = bcell->index();
            facet_computed[id] = true;
            Facet facet(mesh, boundary.facet_index(*bcell));
            Cell cell(mesh, facet.entities(tdim)[0]);
            uint const local_facet = cell.index(facet);
            facet_data[facet_data_size * id] = cell.facet_area(local_facet);
            cell.normal(local_facet , &facet_data[facet_data_size * id + 1]);
            ++facet_count;
          }
        }
        message("Number of selected %u facets.", facet_count);

        // Wait for completion of transfer and cleanup unneeded data
        message("Wait for completion of transfer and cleanup unneeded data");
        MPI_Waitall(adjs.size(), &sendreq_u[0],&status[0]);
        delete[] sendbuf_u;
        delete[] sendsize_u;
        delete[] sendreq_u;
        MPI_Waitall(adjs.size(), &recvreq_u[0],&status[0]);
        delete[] recvreq_u;

        // Construct vertex - facets map from adjacent rank data
        message("Construct vertex - facets map from adjacent rank data");
        tic();
        _map<uint, Array<real *> > vertex_facets;
        for (uint i = 0; i < adjs.size(); ++i)
        {
          uint const a = adjranks[i];
          uint const size = recvsize_u[i];
          //message("adj%8u: size = %u", a, size);
          uint * const ubuffer = &recvbuf_u[offsets_u[i]];
          real * const rbuffer = &recvbuf_r[offsets_r[i]];
          for (uint k = 0; k < size; k += (2 + ubuffer[k + 1]))
          {
            //message("\tv%8u:", ubuffer[k]);
            for (uint f = 1; f <= ubuffer[k + 1]; ++f)
            {
              uint const ii = ubuffer[k + 1 + f] * facet_data_size;
              dolfin_assert(ii < recvsize_r[i]);
              //message("\t\tf%8u:", ubuffer[k + 1 + f]);
              vertex_facets[ubuffer[k]].push_back(rbuffer + ii);
            }
          }
        }
        delete[] recvbuf_u;
        delete[] offsets_u;
        delete[] recvsize_u;
        tocd();

        // Wait for completion of transfer and cleanup unneeded data
        message("Wait for completion of transfer and cleanup unneeded data");
        MPI_Waitall(adjs.size(), &sendreq_r[0],&status[0]);
        delete[] sendbuf_r;
        delete[] sendsize_r;
        delete[] sendreq_r;
        MPI_Waitall(adjs.size(), &recvreq_r[0],&status[0]);
        delete[] recvreq_r;

        //
        message("Compute vertex normals");
        tic();

        for (_map<uint, Array<real *> >::const_iterator v_it = vertex_facets.begin();
             v_it != vertex_facets.end(); ++v_it)
        {
          //message("v%8u", v_it->first);
          real w = 0.0;
          real n[3] = { 0.0 };
          for (Array<real *>::const_iterator f_it = v_it->second.begin();
               f_it != v_it->second.end(); ++f_it)
          {
            real * const data = (*f_it);

            /*
             * Compute normal
             */

            //message("\tw : %+8f @%p", data[0], data);
            for (uint d = 1; d <= gdim; ++d)
            {
            //  message("\tn%u: %+8f", d, data[d]);
              w += data[0];
              n[d] += data[0] * data[d];
            }
            dolfin_assert(w > 0.0);
          }
        }
        tocd();
        delete[] recvbuf_r;
        delete[] offsets_r;
        delete[] recvsize_r;

        //
        delete[] status;
        delete[] adjranks;
      }
      else
      {
        // Compute vertex normal for inner vertices
        for (CellIterator bcell(boundary); !bcell.end(); ++bcell)
        {
          if ((subdomain == NULL) || (subdomain->inside(*bcell, on_boundary)))
          {
            uint const id = bcell->index();
            Facet facet(mesh, boundary.facet_index(*bcell));
            Cell cell(mesh, facet.entities(tdim)[0]);
            uint const local_facet = cell.index(facet);
            facet_data[facet_data_size * id] = cell.facet_area(local_facet);
            cell.normal(local_facet , &facet_data[facet_data_size * id + 1]);
            ++facet_count;
          }
          //
        }
      }

      delete [] facet_data;
      delete [] facet_computed;

    }
    T.end();
  }
  //---------------------------------------------------------------------------
  return 0;
}
