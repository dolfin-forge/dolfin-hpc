// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

// Modified by Martin Alnes, 2008

// First added:  2007-03-01
// Last changed: 2008-04-10

#include <dolfin/common/types.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Vertex.h>
#include "UFCCell.h"
#include "DofMap.h"
#include "SubSystem.h"
#include <dolfin/common/Array.h>
#include <dolfin/elements/ElementLibrary.h>
#include "UFC.h"
#include <dolfin/main/MPI.h>

#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/GlobalFacetMap.h>
#include <string.h>

#ifdef HAS_MPI
#include <mpi.h>
#endif
using namespace dolfin;

//-----------------------------------------------------------------------------
DofMap::DofMap(ufc::dof_map& dof_map, Mesh& mesh, bool dof_map_local) : dof_map(0), 
               ufc_dof_map(&dof_map), ufc_dof_map_local(false), 
               dolfin_mesh(mesh), num_cells(mesh.numCells()), 
               partitions(0)
{
  // Assume responsibilty for ufc_dof_map
  if(dof_map_local) 
    ufc_dof_map_local = ufc_dof_map;
  init();
}
//-----------------------------------------------------------------------------
DofMap::DofMap(ufc::dof_map& dof_map, Mesh& mesh, MeshFunction<uint>& partitions,
               bool dof_map_local) : dof_map(0), ufc_dof_map(&dof_map), 
               ufc_dof_map_local(false), dolfin_mesh(mesh), num_cells(mesh.numCells()), 
               partitions(&partitions)
{
  // Assume responsibilty for ufc_dof_map
  if(dof_map_local) 
    ufc_dof_map_local = ufc_dof_map;
  init();
}
//-----------------------------------------------------------------------------
DofMap::DofMap(const std::string signature, Mesh& mesh) 
  : dof_map(0), ufc_dof_map(0), ufc_dof_map_local(false),
    dolfin_mesh(mesh), num_cells(mesh.numCells()), partitions(0)
{
  // Create ufc dof map from signature
  ufc_dof_map = ElementLibrary::create_dof_map(signature);
  if (!ufc_dof_map)
    error("Unable to find dof map in library: \"%s\".",signature.c_str());

  // Take resposibility for ufc dof map
  ufc_dof_map_local = ufc_dof_map;

  init();
}
//-----------------------------------------------------------------------------
DofMap::DofMap(const std::string signature, Mesh& mesh, 
               MeshFunction<uint>& partitions) 
  : dof_map(0), ufc_dof_map(0), 
    ufc_dof_map_local(false), dolfin_mesh(mesh), num_cells(mesh.numCells()),
    partitions(&partitions)
{
  // Create ufc dof map from signature
  ufc_dof_map = ElementLibrary::create_dof_map(signature);
  if (!ufc_dof_map)
    error("Unable to find dof map in library: \"%s\".",signature.c_str());

  // Take resposibility for ufc dof map
  ufc_dof_map_local = ufc_dof_map;

  init();
}
//-----------------------------------------------------------------------------
DofMap::~DofMap()
{
  if (dof_map)
  {
    for(uint i = 0; i < dolfin_mesh.numCells(); i++)
      delete[] dof_map[i];
    delete [] dof_map;
  }

  if (ufc_dof_map_local)
    delete ufc_dof_map_local;
}
//-----------------------------------------------------------------------------
DofMap* DofMap::extractDofMap(const Array<uint>& sub_system, uint& offset) const
{
  // Check that dof map has not be re-ordered
  //  if (dof_map)
  //    error("Dof map has been re-ordered. Don't yet know how to extract sub dof maps.");

  // Reset offset
  offset = 0;

  // Recursively extract sub dof map
  ufc::dof_map* sub_dof_map = extractDofMap(*ufc_dof_map, offset, sub_system);
  message(2, "Extracted dof map for sub system: %s", sub_dof_map->signature());
  message(2, "Offset for sub system: %d", offset);

  if (partitions)
    return new DofMap(*sub_dof_map, dolfin_mesh, *partitions, true);
  else
    return new DofMap(*sub_dof_map, dolfin_mesh, true);
}
//-----------------------------------------------------------------------------
ufc::dof_map* DofMap::extractDofMap(const ufc::dof_map& dof_map, uint& offset, const Array<uint>& sub_system) const
{
  // Check if there are any sub systems
  if (dof_map.num_sub_dof_maps() == 0)
    error("Unable to extract sub system (there are no sub systems).");

  // Check that a sub system has been specified
  if (sub_system.size() == 0)
    error("Unable to extract sub system (no sub system specified).");
  
  // Check the number of available sub systems
  if (sub_system[0] >= dof_map.num_sub_dof_maps())
    error("Unable to extract sub system %d (only %d sub systems defined).",
                  sub_system[0], dof_map.num_sub_dof_maps());

  // Add to offset if necessary
  for (uint i = 0; i < sub_system[0]; i++)
  {
    ufc::dof_map* ufc_dof_map = dof_map.create_sub_dof_map(i);
    // FIXME: Can we avoid creating a DofMap here just for getting the global dimension?
    if(partitions)
      DofMap dof_map_test(*ufc_dof_map, dolfin_mesh, *partitions);
    else
      DofMap dof_map_test(*ufc_dof_map, dolfin_mesh);
    offset += ufc_dof_map->global_dimension();
    delete ufc_dof_map;
  }
  
  // Create sub system
  ufc::dof_map* sub_dof_map = dof_map.create_sub_dof_map(sub_system[0]);
  
  // Return sub system if sub sub system should not be extracted
  if (sub_system.size() == 1)
    return sub_dof_map;

  // Otherwise, recursively extract the sub sub system
  Array<uint> sub_sub_system;
  for (uint i = 1; i < sub_system.size(); i++)
    sub_sub_system.push_back(sub_system[i]);
  ufc::dof_map* sub_sub_dof_map = extractDofMap(*sub_dof_map, offset, sub_sub_system);
  delete sub_dof_map;

  return sub_sub_dof_map;
}
//-----------------------------------------------------------------------------
void DofMap::init()
{
  //dolfin_debug("Initializing dof map...");

  // Order vertices, so entities will be created correctly according to convention
  dolfin_mesh.order();

  // Initialize mesh entities used by dof map
  for (uint d = 0; d <= dolfin_mesh.topology().dim(); d++)
    if ( ufc_dof_map->needs_mesh_entities(d) )
      dolfin_mesh.init(d);
  
  // Initialize UFC mesh data (must be done after entities are created)
  ufc_mesh.init(dolfin_mesh);

  // Initialize UFC dof map
  const bool init_cells = ufc_dof_map->init_mesh(ufc_mesh);
  if ( init_cells )
  {
    CellIterator cell(dolfin_mesh);
    UFCCell ufc_cell(*cell);
    for (; !cell.end(); ++cell)
    {
      ufc_cell.update(*cell, dolfin_mesh.distdata());
      ufc_dof_map->init_cell(ufc_mesh, ufc_cell);
    }
    ufc_dof_map->init_cell_finalize();
  }

  //dolfin_debug("Dof map initialized");
}
//-----------------------------------------------------------------------------
void DofMap::tabulate_dofs(uint* dofs, ufc::cell& ufc_cell, uint cell_index)
{
  // Either lookup pretabulated values (if build() has been called)
  // or ask the ufc::dof_map to tabulate the values

  if (dof_map)
  {
    /*
    for (uint i = 0; i < local_dimension(); i++)
      dofs[i] = dof_map[cell_index][i];    
    */
    memcpy(dofs, dof_map[cell_index], sizeof(uint)*local_dimension());
    //memcpy(dofs, dof_map[cell_index], sizeof(uint)*local_dimension()); // FIXME: Maybe memcpy() can be used to speed this up? Test this!
  }
  else
    ufc_dof_map->tabulate_dofs(dofs, ufc_mesh, ufc_cell);
} 
//-----------------------------------------------------------------------------
void DofMap::tabulate_dofs(uint* dofs, const ufc::cell& ufc_cell, uint cell_index) const
{
  // Either lookup pretabulated values (if build() has been called)
  // or ask the ufc::dof_map to tabulate the values
  if (dof_map)
  {
    /*
    for (uint i = 0; i < local_dimension(); i++)
      dofs[i] = dof_map[cell_index][i];
    */
    memcpy(dofs, dof_map[cell_index], sizeof(uint)*local_dimension());
    //memcpy(dofs, dof_map[cell_index], sizeof(uint)*local_dimension()); // FIXME: Maybe memcpy() can be used to speed this up? Test this!
  }
  else
    ufc_dof_map->tabulate_dofs(dofs, ufc_mesh, ufc_cell);
} 
//-----------------------------------------------------------------------------
void DofMap::build(UFC& ufc, uint jj)
{

  if( dof_map ) {
    for(uint i = 0; i < dolfin_mesh.numCells(); i++)
      delete[] dof_map[i];
    delete [] dof_map;
  }

  map.clear();
  // delete[] dof_map;
  //    return;

  if(MPI::numProcesses() == 1) {
     uint *dofs =  new uint[local_dimension()];
     
     dof_map = new uint*[dolfin_mesh.numCells()];

     uint num = 0;

     uint tt = local_dimension() / ufc_dof_map->geometric_dimension();
     for(CellIterator c(dolfin_mesh); !c.end(); ++c) {
       
       dof_map[c->index()] = new uint[local_dimension()];    
       ufc.update(*c, dolfin_mesh.distdata());
            
       for(uint j = 0; j < ufc.form.rank(); j++) {
	 ufc_dof_map->tabulate_dofs(dofs, ufc.mesh, ufc.cell);      

	 if(tt == 0) {
	   for(uint i = 0; i < local_dimension(); i++) {
	     const uint dof = dofs[i];

	     std::map<uint, uint>::iterator it = map.find(dof);
	     if(it != map.end())
	       dof_map[c->index()][i] = it->second;
	     else {
	       dof_map[c->index()][i] = num;
	       map[dof] = num++;
#ifdef BLOCKED
	       dof_map[c->index()][i] /= 2;
	       map[dof] /= 2;
#endif
	     }
	   }	     
	 }
	 else 
	   for(uint i = 0; i < tt; i++) {
	     for(uint k = 0; k < ufc_dof_map->geometric_dimension(); k++) {
	       const uint dof = dofs[i + k * tt];
	       
	       std::map<uint, uint>::iterator it = map.find(dof);
	       if (it != map.end())
		 dof_map[c->index()][i + k * tt] = it->second;
	       else {
		 dof_map[c->index()][i + k * tt] = num;
		 map[dof] = num++;
#ifdef BLOCKED
		 dof_map[c->index()][i + k * tt] /= 2;
		 map[dof] /= 2;
#endif
	       }
	     }
	   }    
       }
     }
     

//             for(CellIterator c(dolfin_mesh); !c.end(); ++c) {
//         ufc.update(*c);
//          for(uint j = 0; j< ufc.form.rank(); j++) {
//   	tabulate_dofs(dofs, ufc.cell, c->index());
//       for(uint i = 0; i < local_dimension(); i++)
//       cout<< dofs[i] << " ";
//       cout<<endl;
// 	 }


	 //   }
     delete[] dofs;
  }
  else { 
#ifdef HAS_MPI
    uint local_dim = local_dimension();
    uint *dofs =  new uint[local_dimension()];

    uint pe_size = MPI::numProcesses();
    uint rank = MPI::processNumber();


    if (ufc_dof_map->global_dimension() == dolfin_mesh.distdata().global_numVertices()) {

      // Make sure the mesh is lineary numbered
      dolfin_mesh.renumber();

      dof_map = new uint*[dolfin_mesh.numCells()];      

      for(CellIterator c(dolfin_mesh); !c.end(); ++c) {
	dof_map[c->index()] = new uint[local_dim];    
	uint i = 0;
	for(VertexIterator v(*c); !v.end(); ++v)
	  dof_map[c->index()][i++] = dolfin_mesh.distdata().get_global(*v);
      }
      
      delete[] dofs;
    }
    else if(ufc_dof_map->global_dimension() == 
       ufc_dof_map->geometric_dimension() * dolfin_mesh.distdata().global_numVertices()) {

      // Make sure the mesh is lineary numbered
      dolfin_mesh.renumber();
      
      uint gdim = ufc_dof_map->geometric_dimension();
      uint num_local = dolfin_mesh.numVertices() - dolfin_mesh.distdata().num_ghost(0);
      
      uint num_dofs = gdim * num_local;
      uint offset = 0;

#if ( MPI_VERSION > 1 )
      MPI_Exscan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else 
      MPI_Scan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      offset -= num_dofs;
#endif
      
      _map<uint,uint> v_offset;

      for(VertexIterator v(dolfin_mesh); !v.end(); ++v) {
	if(!dolfin_mesh.distdata().is_ghost(v->index(), 0)) {
	  v_offset[dolfin_mesh.distdata().get_global(*v)] = offset; 
	  offset += gdim;
	}
      }

      Array<uint> *ghost_buff = new Array<uint>[pe_size];
      for(MeshGhostIterator iter(dolfin_mesh.distdata(), 0); !iter.end(); ++iter)
	ghost_buff[iter.owner()].push_back(dolfin_mesh.distdata().get_global(iter.index(), 0)); 
		
      MPI_Status status;
      Array<uint> send_buff;
      uint src,dest;
      uint recv_size = dolfin_mesh.distdata().num_ghost(0); 
      int recv_count, recv_size_gh, send_size;  
      
      for(uint i = 0; i < pe_size; i++) {
	send_size = ghost_buff[i].size();
	MPI_Reduce(&send_size, &recv_size_gh, 1, 
		   MPI_INT, MPI_SUM, i, MPI_COMM_WORLD);
      }
      
      uint *recv_ghost = new uint[ recv_size_gh];
      uint *recv_buff = new uint[ recv_size ];
      
      for(uint j=1; j < pe_size; j++){
	src = (rank - j + pe_size) % pe_size;
	dest = (rank + j) % pe_size;
	
	MPI_Sendrecv(&ghost_buff[dest][0], ghost_buff[dest].size(),
		     MPI_UNSIGNED, dest, 1, recv_ghost, recv_size_gh, 
		     MPI_UNSIGNED, src, 1, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);
	
	for(int k=0; k < recv_count; k++)
	  send_buff.push_back(v_offset[recv_ghost[k]]);
	
	MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, src, 2,
		     recv_buff, recv_size , MPI_UNSIGNED, dest, 2, 
		     MPI_COMM_WORLD,&status);
	MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);

	for(int j=0; j < recv_count; j++)
	  v_offset[ghost_buff[dest][j]] = recv_buff[j];

	send_buff.clear();
      }

      delete[] recv_ghost;
      delete[] recv_buff;

      dof_map = new uint*[dolfin_mesh.numCells()];      

      for(CellIterator c(dolfin_mesh); !c.end(); ++c) {	
	
	dof_map[c->index()] = new uint[local_dim];    

	uint j = 0;
	for(uint i = 0; i < gdim; i++) {
	  for(VertexIterator v(*c); !v.end(); ++v) {
	    dof_map[c->index()][j++] = v_offset[dolfin_mesh.distdata().get_global(*v)] + i;
	  }
	}

      }

      delete[] dofs;

      for(uint i =0; i < pe_size; i++)
	ghost_buff[i].clear();
      delete[] ghost_buff;
    }
    else {

    
      BoundaryMesh local_boundary;
      local_boundary.init_interior(dolfin_mesh);
      
      dolfin_assert(local_boundary.size(0) > 0);
            
      MeshFunction<uint>* cell_map = local_boundary.data().meshFunction("cell map");
      
      Array<uint> send_buff, send_buff_id;
      std::set<uint> shared_dofs, forbidden_dof;
      std::map<uint, uint> dof_vote;
      
    
      for(CellIterator bc(local_boundary); !bc.end(); ++bc) {
	Facet f(dolfin_mesh, cell_map->get(*bc));
	
	for(CellIterator c(f); !c.end(); ++c) {
	  ufc.update(*c, dolfin_mesh.distdata());    
	  
	  for(uint j =0 ; j < ufc.form.rank(); j++) {
	    ufc_dof_map->tabulate_dofs(dofs, ufc.mesh, ufc.cell);      
	    for(uint i = 0; i < local_dim; i++) {
	      const uint dof = dofs[i];
	      
	      if( shared_dofs.count(dof) == 0 ) {
	      forbidden_dof.insert( dof );
	      shared_dofs.insert( dof );
	      dolfin_assert(dof_vote.count(dof) == 0);
	      dof_vote[ dof ] = rank;
	      send_buff.push_back( dof );	
	      send_buff_id.push_back( dof_vote[dof] );
	      }
	    }
	  }
	}
      }
      
      MPI_Status status;
      int recv_count;
      uint src, dest, num_glb, num_sdof;
      num_sdof = send_buff.size();
      MPI_Allreduce(&num_sdof, &num_glb, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      
      uint *recv_buff = new uint[num_glb];
      uint *recv_buff_id = new uint[num_glb];
      
      for(int k = 1 ; k < (int) pe_size; k++){
	
	src = (rank - k + pe_size) % pe_size;
	dest = (rank + k) % pe_size;    
	
	MPI_Sendrecv(&send_buff_id[0], num_sdof , MPI_UNSIGNED, dest, 1, 
		     recv_buff_id, num_glb , MPI_UNSIGNED, src, 1, 
		     MPI_COMM_WORLD, &status);
	
	MPI_Sendrecv(&send_buff[0], num_sdof , MPI_UNSIGNED, dest, 1, 
		     recv_buff, num_glb , MPI_UNSIGNED, src, 1, 
		     MPI_COMM_WORLD, &status);
	MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  
	
	for(int i = 0; i < recv_count; i++) {
	  if( shared_dofs.count(recv_buff[i]) > 0) {
	    dolfin_assert( dof_vote.count(recv_buff[i]) );
	    if( recv_buff_id[i] < dof_vote[recv_buff[i]] ||
		(recv_buff_id[i] == dof_vote[recv_buff[i]] &&
		 status.MPI_SOURCE < (int) rank)) 
	      shared_dofs.erase(recv_buff[i]);
	  }
	}
      }
      
      for(CellIterator c(dolfin_mesh); !c.end(); ++c) {
	
	ufc.update(*c, dolfin_mesh.distdata());
	
	for(uint j = 0; j < ufc.form.rank(); j++) {
	  ufc_dof_map->tabulate_dofs(dofs, ufc.mesh, ufc.cell);      
	  for(uint i = 0; i < local_dim; i++) {  
	    const uint dof = dofs[i];
	    
	  if(forbidden_dof.count(dof)) 
	    continue;
	  
	  shared_dofs.insert( dof );
	  }
	}
      } 
      
      // Initialize range for each processor
      uint offset = 0;
      uint range = shared_dofs.size();

#if ( MPI_VERSION > 1 )
      MPI_Exscan(&range, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else 
      MPI_Scan(&range, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      offset -= range;
#endif    
      map.clear();
      
      _map<uint, Array<std::pair<uint, uint> > > cell_dof;
  
      send_buff.clear();
      send_buff_id.clear();
      
      dof_map = new uint*[dolfin_mesh.numCells()];
      
      for(CellIterator c(dolfin_mesh); !c.end(); ++c) {
	
	dof_map[c->index()] = new uint[local_dim];    
      
	ufc.update(*c, dolfin_mesh.distdata());
	
	for(uint j = 0; j < ufc.form.rank(); j++) {
	  ufc_dof_map->tabulate_dofs(dofs, ufc.mesh, ufc.cell);      
	  for(uint i = 0; i < local_dim; i++) {
	    const uint dof = dofs[i];
	    dof_map[c->index()][i] = 1;        
	    if(forbidden_dof.count(dof) && (shared_dofs.count(dof) == 0)) {
	      std::pair<uint, uint> row_dof(i, dof);
	      cell_dof[c->index()].push_back(row_dof);
	    continue;
	    }
	    
	    std::map<uint, uint>::iterator it = map.find(dof);
	    if (it != map.end()) {
	      dof_map[c->index()][i] = it->second;
	    }
	    else {
	      dof_map[c->index()][i] = offset; 
	      map[dof] = offset++;
	      if( shared_dofs.count(dof) ) {
		send_buff.push_back( dof );
		send_buff_id.push_back( map[dof] );
	      }
	    }     
	  }    
	}    
      }
      delete[] recv_buff_id;
      delete[] recv_buff;
      
      num_sdof = send_buff.size();
      MPI_Allreduce(&num_sdof, &num_glb, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      
      recv_buff = new uint[num_glb];
      recv_buff_id = new uint[num_glb];
      
      for(int j = 1 ; j < (int) pe_size; j++){
	
	src = (rank - j + pe_size) % pe_size;
	dest = (rank + j) % pe_size;    
	
	MPI_Sendrecv(&send_buff_id[0], num_sdof , MPI_UNSIGNED, dest, 1, 
		     recv_buff_id, num_glb , MPI_UNSIGNED, src, 1, 
		     MPI_COMM_WORLD, &status);
	
	MPI_Sendrecv(&send_buff[0], num_sdof , MPI_UNSIGNED, dest, 1, 
		     recv_buff, num_glb , MPI_UNSIGNED, src, 1, 
		   MPI_COMM_WORLD, &status);
	MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);  
	
	for(int i = 0; i < recv_count; i++)  {
	  dolfin_assert( !map.count(recv_buff[i]) );
	  if(shared_dofs.find(recv_buff[i]) != shared_dofs.end())
	    map[ recv_buff[i] ] = recv_buff_id[i];
	}
	
      }
      delete[] recv_buff_id;
      delete[] recv_buff;
      
      
      _map< uint, Array<std::pair<uint, uint> >  >::iterator cit;
      std::vector< std::pair<uint, uint> >::iterator rit;
      for(cit = cell_dof.begin(); cit != cell_dof.end(); ++cit) 
	for(rit = cit->second.begin(); rit != cit->second.end(); ++rit)
	  dof_map[cit->first][rit->first] = map[rit->second];
      
      delete[] dofs;
    }
#endif
  }
}
//-----------------------------------------------------------------------------
std::map<dolfin::uint, dolfin::uint> DofMap::getMap() //FIXME: const
{
  return map;
}
//-----------------------------------------------------------------------------
void DofMap::disp() const
{
  cout << "DofMap" << endl;
  cout << "------" << endl;
  
  // Begin indentation
  begin("");

  // Display UFC dof_map information
  cout << "ufc::dof_map info" << endl;
  cout << "-----------------" << endl;
  begin("");

  cout << "Signature:            " << ufc_dof_map->signature() << endl;
  cout << "Global dimension:     " << ufc_dof_map->global_dimension() << endl;
  cout << "Local dimension:      " << ufc_dof_map->local_dimension() << endl;
  cout << "Geometric dimension:  " << ufc_dof_map->geometric_dimension() << endl;
  cout << "Number of subdofmaps: " << ufc_dof_map->num_sub_dof_maps() << endl;
  cout << "Number of facet dofs: " << ufc_dof_map->num_facet_dofs() << endl;

  for(uint d=0; d<=dolfin_mesh.topology().dim(); d++)
  {
    cout << "Number of entity dofs (dim " << d << "): " << ufc_dof_map->num_entity_dofs(d) << endl;
  }
  for(uint d=0; d<=dolfin_mesh.topology().dim(); d++)
  {
    cout << "Needs mesh entities (dim " << d << "):   " << ufc_dof_map->needs_mesh_entities(d) << endl;
  }
  cout << endl;
  end();

  // Display mesh information
  cout << "Mesh info" << endl;
  cout << "---------" << endl;
  begin("");
  cout << "Geometric dimension:   " << dolfin_mesh.geometry().dim() << endl;
  cout << "Topological dimension: " << dolfin_mesh.topology().dim() << endl;
  cout << "Number of vertices:    " << dolfin_mesh.numVertices() << endl;
  cout << "Number of edges:       " << dolfin_mesh.numEdges() << endl;
  cout << "Number of faces:       " << dolfin_mesh.numFaces() << endl;
  cout << "Number of facets:      " << dolfin_mesh.numFacets() << endl;
  cout << "Number of cells:       " << dolfin_mesh.numCells() << endl;
  cout << endl;
  end();

  cout << "Local cell dofs associated with cell entities (tabulate_entity_dofs output):" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
  begin("");
  {
    uint tdim = dolfin_mesh.topology().dim();
    for(uint d=0; d<=tdim; d++)
    {
      uint num_dofs = ufc_dof_map->num_entity_dofs(d);
      if(num_dofs)
      {
        uint num_entities = dolfin_mesh.type().numEntities(d);
        uint* dofs = new uint[num_dofs];
        for(uint i=0; i<num_entities; i++)
        {
          cout << "Entity (" << d << ", " << i << "):  ";
          ufc_dof_map->tabulate_entity_dofs(dofs, d, i);
          for(uint j=0; j<num_dofs; j++)
          {
            cout << dofs[j];
            if(j < num_dofs-1) cout << ", ";
          }
          cout << endl;
        }
        delete [] dofs;
      }
    }
    cout << endl;
  }
  end();

  cout << "Local cell dofs associated with facets (tabulate_facet_dofs output):" << endl;
  cout << "--------------------------------------------------------------------" << endl;
  begin("");
  {
    uint tdim = dolfin_mesh.topology().dim();
    uint num_dofs = ufc_dof_map->num_facet_dofs();
    uint num_facets = dolfin_mesh.type().numEntities(tdim-1);
    uint* dofs = new uint[num_dofs];
    for(uint i=0; i<num_facets; i++)
    {
      cout << "Facet " << i << ":  ";
      ufc_dof_map->tabulate_facet_dofs(dofs, i);
      for(uint j=0; j<num_dofs; j++)
      {
        cout << dofs[j];
        if(j < num_dofs-1) cout << ", ";
      }
      cout << endl;
    }
    delete [] dofs;
    cout << endl;
  }
  end();

  cout << "tabulate_dofs output" << endl;
  cout << "--------------------" << endl;
  begin("");
  {
    uint tdim = dolfin_mesh.topology().dim();
    uint num_dofs = ufc_dof_map->local_dimension();
    uint* dofs = new uint[num_dofs];
    CellIterator cell(dolfin_mesh);
    UFCCell ufc_cell(*cell);
    for (; !cell.end(); ++cell)
    {
      ufc_cell.update(*cell, dolfin_mesh.distdata());
 
      ufc_dof_map->tabulate_dofs(dofs, ufc_mesh, ufc_cell);
 
      cout << "Cell " << ufc_cell.entity_indices[tdim][0] << ":  ";
      for(uint j=0; j<num_dofs; j++)
      {
        cout << dofs[j];
        if(j < num_dofs-1) cout << ", ";
      }
      cout << endl;
    }
    delete [] dofs;
    cout << endl;
  }
  end();

  cout << "tabulate_coordinates output" << endl;
  cout << "---------------------------" << endl;
  begin("");
  {
    uint tdim = dolfin_mesh.topology().dim();
    uint gdim = ufc_dof_map->geometric_dimension();
    uint num_dofs = ufc_dof_map->local_dimension();
    double** coordinates = new double*[num_dofs];
    for(uint k=0; k<num_dofs; k++)
    {
      coordinates[k] = new double[gdim];
    }
    CellIterator cell(dolfin_mesh);
    UFCCell ufc_cell(*cell);
    for (; !cell.end(); ++cell)
    {
      ufc_cell.update(*cell, dolfin_mesh.distdata());

      ufc_dof_map->tabulate_coordinates(coordinates, ufc_cell);

      cout << "Cell " << ufc_cell.entity_indices[tdim][0] << ":  ";
      for(uint j=0; j<num_dofs; j++)
      {
        cout << "(";
        for(uint k=0; k<gdim; k++)
        {
          cout << coordinates[j][k];
          if(k < gdim-1) cout << ", ";
        }
        cout << ")";
        if(j < num_dofs-1) cout << ",  ";
      }
      cout << endl;
    }
    for(uint k=0; k<gdim; k++)
    {
      delete [] coordinates[k];
    }
    delete [] coordinates;
    cout << endl;
  }
  end();

  // TODO: Display information on renumbering?
  // TODO: Display information on parallel stuff?
  
  // End indentation
  end();
}
//-----------------------------------------------------------------------------

