// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-02-25
// Last changed: 2008-02-25

#include "Mesh.h"
#include "Facet.h"
#include "Vertex.h"
#include "GlobalFacetMap.h"
#include <dolfin/main/MPI.h>

#include <mpi.h>


using namespace dolfin;
//-----------------------------------------------------------------------------
GlobalFacetMap::GlobalFacetMap(Mesh& mesh) : _mesh(mesh)
{
}
//-----------------------------------------------------------------------------
GlobalFacetMap::~GlobalFacetMap()
{
}
//-----------------------------------------------------------------------------
void GlobalFacetMap::init()
{

  // Generate facet - cell connectivity if not generated
  _mesh.init(_mesh.topology().dim() - 1, _mesh.topology().dim());
  
  // Iterate over all Facets connected to the shared vertices
  for(MeshSharedIterator sv(_mesh.distdata(), 0); !sv.end(); ++sv) {
    Vertex v(_mesh, sv.index());
    for(FacetIterator f(v); !f.end(); ++f) {
      if (f->numEntities(_mesh.topology().dim()) == 1) {
      // Mark all facets as local facets
	if(global_facet.count(f->index()) == 0) 
	  global_facet[f->index()] = false;
      }    
    }
  }

  const uint dim = _mesh.topology().dim();
  
  switch(dim)
    {
    case 2:
      findGlobal2D(); break;
    case 3:
      findGlobal3D(); break;
    default:
      error("Couldnt handle local to global map with facet of dim %d", dim);
    }
}
//-----------------------------------------------------------------------------
void GlobalFacetMap::findGlobal2D()
{
  std::map<uint,bool>::iterator iter;
  
  for(iter = global_facet.begin(); iter != global_facet.end(); ++iter){

    Facet f(_mesh, iter->first);
    uint num_shared = 0;

    for (VertexIterator v(f); !v.end(); ++v) 
      if( _mesh.distdata().is_shared(v->index()) )
	num_shared++;  
    
    if ( f.numEntities(_mesh.topology().dim()) == 1 && 
	   num_shared < _mesh.topology().dim() )
      global_facet[f.index()] = true;
  }

}
//-----------------------------------------------------------------------------
void GlobalFacetMap::findGlobal3D()
{

  Array<uint> send_buff;
  const uint dim = _mesh.topology().dim();
  std::map<uint,bool>::iterator iter;
  std::map<uint, uint>::iterator  uiter;
  std::map<uint, uint> unassigned;
  
  uint num_un = 0;

  for(iter = global_facet.begin(); iter != global_facet.end(); ++iter){
    
    Facet f(_mesh, iter->first);
    uint num_shared = 0;
    
    for (VertexIterator v(f); !v.end(); ++v) 
      if( _mesh.distdata().is_shared(v->index()) )
	num_shared++;  
    
    if ( f.numEntities(_mesh.topology().dim()) == 1 && 
	 num_shared < _mesh.topology().dim() )
      global_facet[f.index()] = true;
    else {
      unassigned[num_un++] = f.index();
      for (VertexIterator v(f); !v.end(); ++v) 
	send_buff.push_back(_mesh.distdata().get_global(*v));
      // Make an optimistic guess
      global_facet[f.index()] = true;
    }    
  }

  int sh_count = send_buff.size();
  int res_size = sh_count/dim;
  int recv_size,recv_count;

  MPI_Allreduce(&sh_count, &recv_size, 1, MPI_INT,MPI_MAX, MPI_COMM_WORLD);
  uint *recv_buff = new uint[ recv_size ];
  uint *res_buff = new uint[ res_size ];

  Array<uint> shared_buff;

  MPI_Status status;
  uint pe_size = MPI::numProcesses();
  uint rank = MPI::processNumber();
  uint src, dest;
  uint num_own = 0;

  for(uint j=1; j<pe_size; j++){
    
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;  

    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, dest, 1,
		 recv_buff, recv_size, MPI_UNSIGNED, src, 1,
		 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  

    for(int i = 0; i < recv_count; i += dim){
      uint vi = 0;
      for(uint k = 0; k < dim; k++) 
	if(! _mesh.distdata().have_global(recv_buff[i+k])) 
      	  vi++;
      if(vi) {
	shared_buff.push_back(0);
	continue;
      }
      
      Vertex vertex(_mesh, _mesh.distdata().get_local(recv_buff[i], 0));

      for(FacetIterator f(vertex); !f.end(); ++f) {	
	num_own = 0;
	// Only consider facets conneted to one cell
	if ( f->numEntities(_mesh.topology().dim()) != 1)
	  continue;
	
	for(VertexIterator v(*f); !v.end(); ++v) {
	  for(uint k = 0; k < dim; k++)
	    if(recv_buff[i+k] == _mesh.distdata().get_global(*v))
	      num_own++;
	}
	if(num_own == dim)
	  break;
      }
      if(num_own == dim)
	shared_buff.push_back(1);
      else
	shared_buff.push_back(0);
    }

    MPI_Sendrecv(&shared_buff[0], shared_buff.size(), MPI_UNSIGNED, src, 2,
		 res_buff, res_size, MPI_UNSIGNED, dest, 2,
		 MPI_COMM_WORLD, &status);
    for(int i = 0; i<res_size; i++) 
      if(res_buff[i] == 1)
	global_facet[unassigned[i]] = false;
    shared_buff.clear();
  }

  delete[] recv_buff;
}
//-----------------------------------------------------------------------------
bool GlobalFacetMap::globalFacet(Facet& facet)
{

  const uint index = facet.index();

  // If the facet is in the map, it might be a local facet
  if(global_facet.count(index) > 0 && MPI::numProcesses() > 1)
    return global_facet[index];
  else
    return (facet.numEntities(_mesh.topology().dim()) == 1);


}
//-----------------------------------------------------------------------------
