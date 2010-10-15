// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
//
// First added:  2007-04-24
// Last changed: 2007-07-21

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/log.h>
#include "MeshEntityIterator.h"
#include "Vertex.h"
#include "SubDomain.h"
#include "GlobalFacetMap.h"
#include <dolfin/main/MPI.h>
#include "Facet.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace dolfin;

//-----------------------------------------------------------------------------
SubDomain::SubDomain()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SubDomain::~SubDomain()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
bool SubDomain::inside(const real* x, bool on_boundary) const
{
  error("Unable to determine if point is inside subdomain, function inside() not implemented by user.");
  return false;
}
//-----------------------------------------------------------------------------
void SubDomain::map(const real* x, real* y) const
{
  error("Mapping between subdomains missing for periodic boundary conditions, function map() not implemented by user.");
}
//-----------------------------------------------------------------------------
void SubDomain::mark(MeshFunction<uint>& sub_domains, uint sub_domain) const
{
  message(1, "Computing sub domain markers for sub domain %d.", sub_domain);

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
  facetmap.init();

  // Always false when not marking facets
  bool on_boundary = false;

  // Compute sub domain markers
  for (MeshEntityIterator entity(mesh, dim); !entity.end(); ++entity)
  {
    on_boundary = false;
    // Check if entity is on the boundary if entity is a facet
    if (dim == D - 1) 
      on_boundary = entity->numEntities(D) == 1;
    else if ( dim == 0)
      for( FacetIterator fi(*entity); !fi.end(); ++fi) {
	/*
	if( fi->numEntities(D) == 1 ||
	    (MPI::numProcesses() > 1 && facetmap.globalFacet(*fi) &&
	     fi->numEntities(D) == 1))
	*/
	if(fi->numEntities(D) &&  facetmap.globalFacet(*fi))
	  on_boundary = true;
      }
    
    bool all_vertices_inside = true;
    // Dimension of facet > 0, check incident vertices
    if (entity->dim() > 0)
    {

      if(MPI::numProcesses() > 1){
	Facet f(mesh, entity->index());
	if(!facetmap.globalFacet(f))
	  on_boundary = false; 
      }

      for (VertexIterator vertex(*entity); !vertex.end(); ++vertex)
      {
        simple_array<real> x(mesh.geometry().dim(), vertex->x());
        if (!inside(x, on_boundary))
        {
          all_vertices_inside = false;
          break;
        }
      }
    }
    // Dimension of facet == 0, so just check the vertex itself
    else
    {

      simple_array<real> x(mesh.geometry().dim(), mesh.geometry().x(entity->index()));
      if (!inside(x, on_boundary))
      {
        all_vertices_inside = false;	
      }
    }

    // Mark entity with all vertices inside
    if (all_vertices_inside)
      sub_domains(*entity) = sub_domain;
  }

#ifdef HAVE_MPI
  if(MPI::numProcesses() > 1) {
    if(dim == 0) {
      uint pe_size = MPI::numProcesses();
      uint rank = MPI::processNumber();
      
      Array<uint> *ghost_buff = new Array<uint>[pe_size];
      for(MeshGhostIterator iter(mesh.distdata(), 0); !iter.end(); ++iter)
	if(sub_domains.get(iter.index()) == sub_domain)
	  ghost_buff[iter.owner()].push_back(mesh.distdata().get_global(iter.index(),0));
      
      int recv_size,recv_count, send_size;
      for(uint i=0; i<pe_size; i++){
	send_size = ghost_buff[i].size();
	MPI_Reduce(&send_size, &recv_size, 1, MPI_INT, MPI_SUM, i, MPI::DOLFIN_COMM);
      }
      uint *recv_buff = new uint[ recv_size];
      
      MPI_Status status;
      uint src,dest;
      
      for(uint j=1; j<pe_size; j++){
	
	src = (rank - j + pe_size) % pe_size;
	dest = (rank + j) % pe_size;    
	
	MPI_Sendrecv(&ghost_buff[dest][0], ghost_buff[dest].size(), MPI_UNSIGNED,
		     dest, 1, recv_buff, recv_size, MPI_UNSIGNED, src, 1,
		     MPI::DOLFIN_COMM, &status);
	MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);      
	
	for(int i = 0; i < recv_count; i++) 
	  sub_domains.set(mesh.distdata().get_local(recv_buff[i], 0), sub_domain);
	
      }
      delete[] recv_buff;
      
      for(uint i = 0; i < pe_size; i++)
	ghost_buff[i].clear();
      delete[] ghost_buff;
      
    }
  }
#endif

}
//-----------------------------------------------------------------------------
