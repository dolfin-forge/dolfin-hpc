// Copyright (C) 2007-2008 Anders Logg and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

// Modified by Martin Alnes, 2008
// Modified by Niclas Jansson, 2009

// First added:  2007-03-01
// Last changed: 2009-11-01

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/common/Array.h>
#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/main/MPI.h>

#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/GlobalFacetMap.h>
#include <cstring>
#include <cstdlib>

#ifdef HAVE_MPI
#include <mpi.h>
#endif
using namespace dolfin;

//-----------------------------------------------------------------------------
DofMap::DofMap(ufc::dof_map& dof_map, Mesh& mesh, bool dof_map_local) : dof_map(0), 
               ufc_dof_map(&dof_map), ufc_dof_map_local(0), 
               dolfin_mesh(mesh), num_cells(mesh.numCells()), 
               partitions(0), _type_(-1), _local_size(0), v_map(0)
{


  // Assume responsibilty for ufc_dof_map
  if(dof_map_local) 
    ufc_dof_map_local = ufc_dof_map;
  init();  

  if(dolfin::MPI::numProcesses() > 1)
    build();
}
//-----------------------------------------------------------------------------
DofMap::DofMap(ufc::dof_map& dof_map, Mesh& mesh, MeshFunction<uint>& partitions,
               bool dof_map_local) : dof_map(0), ufc_dof_map(&dof_map), 
               ufc_dof_map_local(0), dolfin_mesh(mesh), num_cells(mesh.numCells()), 
               partitions(&partitions),  _type_(-1), _local_size(0), v_map(0)
{
  // Assume responsibilty for ufc_dof_map
  if(dof_map_local) 
    ufc_dof_map_local = ufc_dof_map;
  init();

  if(dolfin::MPI::numProcesses() > 1)
    build();
}
//-----------------------------------------------------------------------------
DofMap::DofMap(const std::string signature, Mesh& mesh) 
  : dof_map(0), ufc_dof_map(0), ufc_dof_map_local(0),
    dolfin_mesh(mesh), num_cells(mesh.numCells()), partitions(0), 
    _type_(-1), _local_size(0), v_map(0)
{
  // Create ufc dof map from signature
  ufc_dof_map = ElementLibrary::create_dof_map(signature);
  if (!ufc_dof_map)
    error("Unable to find dof map in library: \"%s\".",signature.c_str());

  // Take resposibility for ufc dof map
  ufc_dof_map_local = ufc_dof_map;

  init();

  if(dolfin::MPI::numProcesses() > 1)
    build();
}
//-----------------------------------------------------------------------------
DofMap::DofMap(const std::string signature, Mesh& mesh, 
               MeshFunction<uint>& partitions) 
  : dof_map(0), ufc_dof_map(0), 
    ufc_dof_map_local(0), dolfin_mesh(mesh), num_cells(mesh.numCells()),
    partitions(&partitions), _type_(-1), _local_size(0), v_map(0)
{
  // Create ufc dof map from signature
  ufc_dof_map = ElementLibrary::create_dof_map(signature);
  if (!ufc_dof_map)
    error("Unable to find dof map in library: \"%s\".",signature.c_str());

  // Take resposibility for ufc dof map
  ufc_dof_map_local = ufc_dof_map;

  init();

  if(dolfin::MPI::numProcesses() > 1)
    build();
}
//-----------------------------------------------------------------------------
DofMap::~DofMap()
{
  if (dof_map)
    delete [] dof_map;

  if (ufc_dof_map_local)
    delete ufc_dof_map_local;

  if(v_map)
    delete[] v_map;
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

  if(_type_ == 0)
  {
    Cell c(dolfin_mesh, cell_index);
    for (uint i = 0; i < local_dimension(); i++)
      dofs[i] = dolfin_mesh.distdata().get_global(c.entities(0)[i], 0);
  }
  else if(_type_ == 1)
  {
    *dofs = dolfin_mesh.distdata().get_cell_global(cell_index);
  }
  else if(_type_ == 2 && v_map)
  {
    Cell c(dolfin_mesh, cell_index);
    uint gdim = ufc_dof_map->num_sub_dof_maps();
    uint num_entities = c.numEntities(0);
    dolfin_assert(gdim * num_entities == local_dimension());
    for (uint k = 0; k < gdim; k++)
      for (uint i = 0;  i < num_entities; i++)
       {
	 dofs[i + k * (num_entities)] = v_map[c.entities(0)[i]] + k;
       }
  }
  else if(_type_ == 3)
  {
    uint gdim = ufc_dof_map->num_sub_dof_maps();
    dolfin_assert(local_dimension() == gdim);
    for(uint i = 0; i < gdim; i++)
      dofs[i] = (cell_index + i * dolfin_mesh.numCells()) + _offset_;
  }
  else if (dof_map)
  {
    uint offset = local_dimension() * cell_index;
    memcpy(dofs, &dof_map[offset], sizeof(uint)*local_dimension());
  }
  else
    ufc_dof_map->tabulate_dofs(dofs, ufc_mesh, ufc_cell);
} 
//-----------------------------------------------------------------------------
void DofMap::tabulate_dofs(uint* dofs, const ufc::cell& ufc_cell, uint cell_index) const
{
  // Either lookup pretabulated values (if build() has been called)
  // or ask the ufc::dof_map to tabulate the values
  if(_type_ == 0)
  {
    Cell c(dolfin_mesh, cell_index);
    for (uint i = 0; i < local_dimension(); i++)
      dofs[i] = dolfin_mesh.distdata().get_global(c.entities(0)[i], 0);
  }
  else if(_type_ == 1)
  {
    *dofs = dolfin_mesh.distdata().get_cell_global(cell_index);
  }
  else if(_type_ == 2 && v_map)
  {
    Cell c(dolfin_mesh, cell_index);
    uint gdim = ufc_dof_map->num_sub_dof_maps();
    uint num_entities = c.numEntities(0);
    for (uint k = 0; k < gdim; k++)
      for (uint i = 0;  i < num_entities; i++)
       {
	 dofs[i + k * (num_entities)] = v_map[c.entities(0)[i]] + k;
       }
  }
  else if(_type_ == 3)
  {
    uint gdim = ufc_dof_map->num_sub_dof_maps();
    for(uint i = 0; i < gdim; i++)
      dofs[i] = (cell_index + i * dolfin_mesh.numCells()) + _offset_;

  }
  else if (dof_map)
  {
    uint offset = local_dimension() * cell_index;
    memcpy(dofs, &dof_map[offset], sizeof(uint)*local_dimension());
  }
  else
    ufc_dof_map->tabulate_dofs(dofs, ufc_mesh, ufc_cell);

}//-----------------------------------------------------------------------------
void DofMap::build()
{

  if( dof_map ) 
    delete [] dof_map;

  map.clear();


  if(MPI::numProcesses() == 1) {
    /*
     uint *dofs =  new uint[local_dimension()];
     
     dof_map = new uint*[dolfin_mesh.numCells()];

     uint num = 0;

     uint tt = local_dimension() / ufc_dof_map->num_sub_dof_maps();
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
	     for(uint k = 0; k < ufc_dof_map->num_sub_dof_maps(); k++) {
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
	 */
    //     delete[] dofs;
  }
  else { 
#ifdef HAVE_MPI
    uint *dofs =  new uint[local_dimension()];

    uint pe_size = MPI::numProcesses();
    uint rank = MPI::processNumber();

    dolfin_mesh.renumber();

    if (ufc_dof_map->global_dimension() == dolfin_mesh.distdata().global_numVertices()) {
      _type_ = 0;
      _local_size = dolfin_mesh.numVertices() - dolfin_mesh.distdata().num_ghost(0);
    }
    else if(ufc_dof_map->global_dimension() == dolfin_mesh.distdata().global_numCells()) {
      _type_ = 1;      
      _local_size = dolfin_mesh.numCells();
    }
    else if(ufc_dof_map->global_dimension() == 
       ufc_dof_map->num_sub_dof_maps() * dolfin_mesh.distdata().global_numVertices()) {
      
      uint gdim = ufc_dof_map->num_sub_dof_maps();
      uint num_local = dolfin_mesh.numVertices() - dolfin_mesh.distdata().num_ghost(0);
      
      uint num_dofs = gdim * num_local;
      uint offset = 0;

#if ( MPI_VERSION > 1 )
      MPI_Exscan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else 
      MPI_Scan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
      offset -= num_dofs;
#endif
      _map<uint, uint> v_offset;

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
                   MPI_INT, MPI_SUM, i, MPI::DOLFIN_COMM);
      }
      
      uint *recv_ghost = new uint[ recv_size_gh];
      uint *recv_buff = new uint[ recv_size ];
      
      for(uint j=1; j < pe_size; j++){
        src = (rank - j + pe_size) % pe_size;
        dest = (rank + j) % pe_size;
        
        MPI_Sendrecv(&ghost_buff[dest][0], ghost_buff[dest].size(),
                     MPI_UNSIGNED, dest, 1, recv_ghost, recv_size_gh, 
                     MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
        MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);
        
        for(int k=0; k < recv_count; k++)
          send_buff.push_back(v_offset[recv_ghost[k]]);
        
        MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, src, 2,
                     recv_buff, recv_size , MPI_UNSIGNED, dest, 2, 
                     MPI::DOLFIN_COMM,&status);
        MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);

        for(int j=0; j < recv_count; j++)
          v_offset[ghost_buff[dest][j]] = recv_buff[j];

        send_buff.clear();
      }

      delete[] recv_ghost;
      delete[] recv_buff;

      if( v_map )
	delete[] v_map;

      v_map = new uint[dolfin_mesh.numVertices()];
      for(VertexIterator v(dolfin_mesh); !v.end(); ++v)
	v_map[v->index()] = v_offset[dolfin_mesh.distdata().get_global(v->index(), 0)];
      
      _type_ = 2;
      v_offset.clear();
      
      _local_size = gdim * num_local;

      for(uint i =0; i < pe_size; i++)
        ghost_buff[i].clear();
      delete[] ghost_buff;
      
    }
    else if(ufc_dof_map->global_dimension() == 
       ufc_dof_map->num_sub_dof_maps() * dolfin_mesh.distdata().global_numCells()) {

      uint gdim = ufc_dof_map->num_sub_dof_maps();
      uint num_dofs =  gdim * dolfin_mesh.numCells();      
      uint offset = 0;

#if ( MPI_VERSION > 1 )
      MPI_Exscan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else 
      MPI_Scan(&num_dofs, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
      offset -= num_dofs;
#endif
      _offset_ = offset;
      _type_ = 3;
      _local_size =  gdim * dolfin_mesh.numCells();      
    }
    else {


      BoundaryMesh interior_boundary;
      interior_boundary.init_interior(dolfin_mesh);
      MeshFunction<uint>* cell_map =interior_boundary.data().meshFunction("cell map");

      std::vector<uint> send_buffer;                                                
      _set<uint> shared_dofs, forbidden_dofs, owned_dofs;                            
      _map<uint, uint> dof_vote;                                                
      _map<uint, std::vector<uint> > dof2index;

      uint n = local_dimension(); 
      dof_map = new uint[n * dolfin_mesh.numCells()];
      uint *facet_dofs = new uint[num_facet_dofs()];

      Cell c_tmp(dolfin_mesh, 1);
      UFCCell ufc_cell(c_tmp);      
                                                                                      
      // Initialize random number generator differently on each process             
      srand((uint)time(0) + MPI::processNumber());   

      // Decide ownership of shared dofs                                            
      for (CellIterator bc(interior_boundary); !bc.end(); ++bc)                     
      {                                                                             
	Facet f(dolfin_mesh, cell_map->get(*bc));
	Cell c(dolfin_mesh, f.entities(dolfin_mesh.topology().dim())[0]);

	uint local_facet = c.index(f);                                          

	ufc_cell.update(c, dolfin_mesh.distdata());                                  
	ufc_dof_map->tabulate_dofs(dofs, ufc_mesh, ufc_cell);                        
	ufc_dof_map->tabulate_facet_dofs(facet_dofs, local_facet);

	for (uint i = 0; i < num_facet_dofs(); i++)    
	{                                                                         
	  // Assign an ownership vote for each "shared" dof                       
	  if (shared_dofs.find(dofs[i]) == shared_dofs.end())                     
	  {                                                                       
	    shared_dofs.insert(dofs[i]);                                          
	    dof_vote[dofs[i]] = (uint) rand() + (uint) MPI::processNumber();
	    send_buffer.push_back(dofs[i]);                                       
	    send_buffer.push_back(dof_vote[dofs[i]]);                             
	  }                                                                       
	}                                                                         
      }  

      // Decide ownership of "shared" dofs                                          
      MPI_Status status;
      int recv_count;
      uint src, dest, max_recv;
      uint num_proc = MPI::numProcesses();                                         
      uint proc_num = MPI::processNumber();                                        
      uint local_size = send_buffer.size();
      MPI_Allreduce(&local_size, &max_recv, 1,
		   MPI_UNSIGNED, MPI_MAX, MPI::DOLFIN_COMM);
      uint *recv_buffer = new uint[max_recv];                                       
      for(uint k = 1; k < MPI::numProcesses(); ++k)                                
      {                                                                             
	src = (proc_num - k + num_proc) % num_proc;                                 
	dest = (proc_num +k) % num_proc;                                            
	
	MPI_Sendrecv(&send_buffer[0], send_buffer.size(), MPI_UNSIGNED, dest, 1,
		     recv_buffer, max_recv, MPI_UNSIGNED, src, 1,
		     MPI::DOLFIN_COMM, &status);
	MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
	
	for (int i = 0; i < recv_count; i += 2)                                    
	{                                                                           
	  if (shared_dofs.find(recv_buffer[i]) != shared_dofs.end())                
	  {                                                                         
	    // Move dofs with higher ownership votes from shared to forbidden       
	    if (recv_buffer[i+1] < dof_vote[recv_buffer[i]] ) {                     
	      forbidden_dofs.insert(recv_buffer[i]);                                
	      shared_dofs.erase(recv_buffer[i]);                                    
	    }                                                                       
	  }                                                                         
	}                                                                           
      }                                                                             
      
      send_buffer.clear();  
      
      // Mark all non forbidden dofs as owned by the processes                      
      for (CellIterator c(dolfin_mesh); !c.end(); ++c)                                     
      {                                                                             
	ufc_cell.update(*c, dolfin_mesh.distdata());                                  
	ufc_dof_map->tabulate_dofs(dofs, ufc_mesh, ufc_cell);  
	for (uint i = 0; i < n; i++)                                                
	{                                                                           
	  if (forbidden_dofs.find(dofs[i]) == forbidden_dofs.end())                 
	  {                                                                         
	    // Mark dof as owned                                                    
	    owned_dofs.insert(dofs[i]);                                             
	  }                                                                         
	  
	  // Create mapping from dof to dof_map offset                              
	  dof2index[dofs[i]].push_back(c->index() * n + i);                         
	}                                                                           
      }                                                                             
      
      // Compute offset for owned and non shared dofs                               
      uint range = owned_dofs.size();                                               
      _local_size = range;
      uint offset = 0;
#if ( MPI_VERSION > 1 )
      MPI_Exscan(&range, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else 
      MPI_Scan(&range, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
      offset -= range;
#endif    

      // Compute renumbering for local and owned shared dofs                        
      for (_set<uint>::iterator it = owned_dofs.begin();                        
	   it != owned_dofs.end(); ++it, offset++)                                  
      {                                                                             
	for(std::vector<uint>::iterator di = dof2index[*it].begin();                
	    di != dof2index[*it].end(); ++di)                                       
	  dof_map[*di] = offset;                                                   
	
	if (shared_dofs.find(*it) != shared_dofs.end())                             
	{                                                                           
	  send_buffer.push_back(*it);                                               
	  send_buffer.push_back(offset);                                            
	}                                                                           
      }                                                                             
      

      // Exchange new dof numbers for shared dofs                                   
      delete[] recv_buffer;                                                         
      local_size = send_buffer.size();
      MPI_Allreduce(&local_size, &max_recv, 1,
		   MPI_UNSIGNED, MPI_MAX, MPI::DOLFIN_COMM);     
      recv_buffer = new uint[max_recv];                                             
      for(uint k = 1; k < MPI::numProcesses(); ++k)                                
      {                                                                             
	src = (proc_num - k + num_proc) % num_proc;                                 
	dest = (proc_num +k) % num_proc;                                            
	
	MPI_Sendrecv(&send_buffer[0], send_buffer.size(), MPI_UNSIGNED, dest, 1,
		     recv_buffer, max_recv, MPI_UNSIGNED, src, 1,
		     MPI::DOLFIN_COMM, &status);
	MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
	
	for (int i = 0; i < recv_count; i += 2)                                    
	{                                                                           
	  // Assign new dof number for shared dofs                                  
	  if (forbidden_dofs.find(recv_buffer[i]) != forbidden_dofs.end())          
	  {                                                                         
	    for(std::vector<uint>::iterator di = dof2index[recv_buffer[i]].begin(); 
		di != dof2index[recv_buffer[i]].end(); ++di)                        
	      dof_map[*di] = recv_buffer[i+1];                                       
	  }                                                                         
	}                                                                           
      }                                                                             
      delete[] recv_buffer;                                                         
      map.clear();
      delete[] facet_dofs;
    }
    delete[] dofs;   

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

  // for(uint d=0; d<=dolfin_mesh.topology().dim(); d++)
  // {
  //   cout << "Number of entity dofs (dim " << d << "): " << ufc_dof_map->num_entity_dofs(d) << endl;
  // }
  // for(uint d=0; d<=dolfin_mesh.topology().dim(); d++)
  // {
  //   cout << "Needs mesh entities (dim " << d << "):   " << ufc_dof_map->needs_mesh_entities(d) << endl;
  // }
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

  // cout << "Local cell dofs associated with cell entities (tabulate_entity_dofs output):" << endl;
  // cout << "----------------------------------------------------------------------------" << endl;
  // begin("");
  // {
  //   uint tdim = dolfin_mesh.topology().dim();
  //   for(uint d=0; d<=tdim; d++)
  //   {
  //     uint num_dofs = ufc_dof_map->num_entity_dofs(d);
  //     if(num_dofs)
  //     {
  //       uint num_entities = dolfin_mesh.type().numEntities(d);
  //       uint* dofs = new uint[num_dofs];
  //       for(uint i=0; i<num_entities; i++)
  //       {
  //         cout << "Entity (" << d << ", " << i << "):  ";
  //         ufc_dof_map->tabulate_entity_dofs(dofs, d, i);
  //         for(uint j=0; j<num_dofs; j++)
  //         {
  //           cout << dofs[j];
  //           if(j < num_dofs-1) cout << ", ";
  //         }
  //         cout << endl;
  //       }
  //       delete [] dofs;
  //     }
  //   }
  //   cout << endl;
  // }
  // end();

  // cout << "Local cell dofs associated with facets (tabulate_facet_dofs output):" << endl;
  // cout << "--------------------------------------------------------------------" << endl;
  // begin("");
  // {
  //   uint tdim = dolfin_mesh.topology().dim();
  //   uint num_dofs = ufc_dof_map->num_facet_dofs();
  //   uint num_facets = dolfin_mesh.type().numEntities(tdim-1);
  //   uint* dofs = new uint[num_dofs];
  //   for(uint i=0; i<num_facets; i++)
  //   {
  //     cout << "Facet " << i << ":  ";
  //     ufc_dof_map->tabulate_facet_dofs(dofs, i);
  //     for(uint j=0; j<num_dofs; j++)
  //     {
  //       cout << dofs[j];
  //       if(j < num_dofs-1) cout << ", ";
  //     }
  //     cout << endl;
  //   }
  //   delete [] dofs;
  //   cout << endl;
  // }
  // end();

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

