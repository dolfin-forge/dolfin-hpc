// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
//
// First added:  2007-04-24
// Last changed: 2007-07-21

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/mesh/GlobalFacetMap.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/parameter/parameters.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace dolfin;

//-----------------------------------------------------------------------------
SubDomain::SubDomain() :
    intersection_detector(NULL),
    BMARG(1.0e-6)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SubDomain::SubDomain(Mesh& mesh) :
    intersection_detector(new IntersectionDetector(mesh)),
    BMARG(1.0e-6)
{
}
//-----------------------------------------------------------------------------
SubDomain::~SubDomain()
{
  delete intersection_detector;
}
//-----------------------------------------------------------------------------
bool SubDomain::inside(real const * x, bool const on_boundary) const
{
  error("Unable to determine if point is inside subdomain, function inside() "
        "not implemented by user.");
  return false;
}
//-----------------------------------------------------------------------------
bool SubDomain::inside(MeshEntity& entity, bool const on_boundary) const
{
  uint ret = entity.num_entities(0);
  for (VertexIterator v(entity); !v.end(); ++v)
  {
    if (this->inside(v->x(), on_boundary))
    {
      --ret;
    }
  }
  return (ret == 0);
}
//-----------------------------------------------------------------------------
bool SubDomain::close(real x, real const xref) const
{
  return (std::fabs(x - xref) < BMARG);
}
//-----------------------------------------------------------------------------
bool SubDomain::intersect(real const * x, uint const dim,
                          bool const on_boundary) const
{
  Point p;
  for (uint i = 0; i < dim; i++)
    p[i] = x[i];

  return intersect(p, on_boundary);
}
//-----------------------------------------------------------------------------
bool SubDomain::intersect(Point const& p, bool const on_boundary) const
{

  Array<uint> cells;
  intersection_detector->overlap(p, cells);
  if (dolfin_get("SubDomain Intersect Boundary"))
    return (cells.size() > 0) && on_boundary;
  else
    return (cells.size() > 0);
}
//-----------------------------------------------------------------------------
void SubDomain::mark(MeshFunction<uint>& sub_domains, uint sub_domain) const
{
  message(1, "Computing sub domain markers for sub domain %d.", sub_domain);

  // Save GTS tolerances
  real gts_tol = dolfin_get("GTS Tolerance");
  real geom_tri_tol = dolfin_get("Geometrical Tolerance Triangle");
  real geom_tet_tol = dolfin_get("Geometrical Tolerance Tetrahedron");

  if (intersection_detector)
  {
    dolfin_set("GTS Tolerance", 1e-6);
    dolfin_set("Geometrical Tolerance Triangle",
               dolfin_get("SubDomain Geometrical Tolerance"));
    dolfin_set("Geometrical Tolerance Tetrahedron",
               dolfin_get("SubDomain Geometrical Tolerance"));
  }

  // Get the dimension of the entities we are marking
  const uint dim = sub_domains.dim();

  // Compute facet - cell connectivity if necessary
  Mesh& mesh = sub_domains.mesh();
  const uint D = mesh.topology().dim();
  if (dim == D - 1)
  {
    mesh.init(D - 1);
    mesh.init(D - 1, D);
  }

  GlobalFacetMap facetmap(mesh);

  // Always false when not marking facets
  bool on_boundary = false;

  // Compute sub domain markers
  bool const is_distributed = mesh.is_distributed();
  for (MeshEntityIterator entity(mesh, dim); !entity.end(); ++entity)
  {
    on_boundary = false;
    // Check if entity is on the boundary if entity is a facet
    if (dim == D - 1)
    {
      on_boundary = entity->num_entities(D) == 1;
    }
    else if (dim == 0)
    {
      for (FacetIterator fi(*entity); !fi.end(); ++fi)
      {
        if (fi->num_entities(D) && facetmap.is_global(*fi))
        {
          on_boundary = true;
        }
      }
    }

    bool all_vertices_inside = true;
    // Dimension of facet > 0, check incident vertices
    if (entity->dim() > 0)
    {

      if (is_distributed)
      {
        Facet f(mesh, entity->index());
        if (!facetmap.is_global(f))
        {
          on_boundary = false;
        }
      }

      for (VertexIterator vertex(*entity); !vertex.end(); ++vertex)
      {
        if (intersection_detector)
        {
          if (!intersect(vertex->point(), on_boundary))
          {
            all_vertices_inside = false;
            break;
          }
        }
        else if (!inside(vertex->x(), on_boundary))
        {
          all_vertices_inside = false;
          break;
        }
      }
    }
    // Dimension of facet == 0, so just check the vertex itself
    else
    {
      real * x = mesh.geometry().x(entity->index());
      if (intersection_detector)
      {
        if (!intersect(x, mesh.geometry().dim(), on_boundary))
          all_vertices_inside = false;
      }
      else if (!inside(x, on_boundary))
      {
        all_vertices_inside = false;
      }
    }

    // Mark entity with all vertices inside
    if (all_vertices_inside)
      sub_domains(*entity) = sub_domain;
  }

  // Reset GTS tolerances
  dolfin_set("GTS Tolerance", gts_tol);
  dolfin_set("Geometrical Tolerance Triangle", geom_tri_tol);
  dolfin_set("Geometrical Tolerance Tetrahedron", geom_tet_tol);

#ifdef HAVE_MPI
  if (is_distributed)
  {
    if (dim == 0)
    {
      uint pe_size = MPI::numProcesses();
      uint rank = MPI::processNumber();

      Array<uint> *ghost_buff = new Array<uint> [pe_size];
      for (GhostIterator iter(mesh.distdata()[0]); !iter.end(); ++iter)
        if (sub_domains.get(iter.index()) == sub_domain)
          ghost_buff[iter.owner()].push_back(
              mesh.distdata()[0].get_global(iter.index()));

      int recv_size, recv_count, send_size;
      for (uint i = 0; i < pe_size; i++)
      {
        send_size = ghost_buff[i].size();
        MPI_Reduce(&send_size, &recv_size, 1, MPI_INT, MPI_SUM, i,
                   MPI::DOLFIN_COMM);
      }
      uint *recv_buff = new uint[recv_size];

      MPI_Status status;
      uint src, dest;

      for (uint j = 1; j < pe_size; j++)
      {

        src = (rank - j + pe_size) % pe_size;
        dest = (rank + j) % pe_size;

        MPI_Sendrecv(&ghost_buff[dest][0], ghost_buff[dest].size(),
                     MPI_UNSIGNED, dest, 1, recv_buff, recv_size, MPI_UNSIGNED,
                     src, 1, MPI::DOLFIN_COMM, &status);
        MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

        for (int i = 0; i < recv_count; i++)
          sub_domains.set(mesh.distdata()[0].get_local(recv_buff[i]),
                          sub_domain);

      }
      delete[] recv_buff;

      for (uint i = 0; i < pe_size; i++)
        ghost_buff[i].clear();
      delete[] ghost_buff;

    }
  }
#endif

}
//-----------------------------------------------------------------------------
