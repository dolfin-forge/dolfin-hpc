// Copyright (C) 2003-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008.
//
// First added:  2003-10-21
// Last changed: 2008-07-01

#include <cstring>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>
#include "PXMLMesh.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace dolfin;

//-----------------------------------------------------------------------------
PXMLMesh::PXMLMesh(Mesh& mesh) : XMLObject(), _mesh(mesh), state(OUTSIDE), f(0), a(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PXMLMesh::~PXMLMesh()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
#ifdef HAVE_MPI
void PXMLMesh::startElement(const xmlChar *name, const xmlChar **attrs)
{
  switch ( state )
  {
  case OUTSIDE:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "mesh") == 0 )
    {
      readMesh(name, attrs);
      state = INSIDE_MESH;
    }
    
    break;

  case INSIDE_MESH:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "vertices") == 0 )
    {
      readVertices(name, attrs);
      state = INSIDE_VERTICES;
    }
    else if ( xmlStrcasecmp(name, (xmlChar *) "cells") == 0 )
    {
      readCells(name, attrs);
      state = INSIDE_CELLS;
    }
    else if ( xmlStrcasecmp(name, (xmlChar *) "data") == 0 )
    {
      state = INSIDE_DATA;
    }

    break;
    
  case INSIDE_VERTICES:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "vertex") == 0 )
      readVertex(name, attrs);

    break;
    
  case INSIDE_CELLS:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "interval") == 0 )
      readInterval(name, attrs);
    else if ( xmlStrcasecmp(name, (xmlChar *) "triangle") == 0 )
      readTriangle(name, attrs);
    else if ( xmlStrcasecmp(name, (xmlChar *) "tetrahedron") == 0 )
      readTetrahedron(name, attrs);
    
    break;

  case INSIDE_DATA:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "meshfunction") == 0 )
    {
      error("Parsing mesh function is not supported in parallel");
      readMeshFunction(name, attrs);
      state = INSIDE_MESH_FUNCTION;
    }
    if ( xmlStrcasecmp(name, (xmlChar *) "array") == 0 )
    {
      error("Parsing array is not supported in parallel");
      readArray(name, attrs);
      state = INSIDE_ARRAY;
    }

    break;

  case INSIDE_MESH_FUNCTION:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "entity") == 0 )
      readMeshEntity(name, attrs);

    break;

  case INSIDE_ARRAY:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "element") == 0 )
      readArrayElement(name, attrs);

    break;

  default:
    ;
  }
}
//-----------------------------------------------------------------------------
void PXMLMesh::endElement(const xmlChar *name)
{
  switch ( state )
  {
  case INSIDE_MESH:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "mesh") == 0 )
    {
      closeMesh();
      state = DONE;
    }
    
    break;
    
  case INSIDE_VERTICES:
    
    if ( xmlStrcasecmp(name, (xmlChar *) "vertices") == 0 )
    {
      state = INSIDE_MESH;    
    }

    break;

  case INSIDE_CELLS:
	 
    if ( xmlStrcasecmp(name, (xmlChar *) "cells") == 0 )
    {
      state = INSIDE_MESH;
    }

    break;

  case INSIDE_DATA:

    if ( xmlStrcasecmp(name, (xmlChar *) "data") == 0 )
    {
      state = INSIDE_MESH;
    }

    break;

  case INSIDE_MESH_FUNCTION:

    if ( xmlStrcasecmp(name, (xmlChar *) "meshfunction") == 0 )
    {
      state = INSIDE_DATA;
    }

    break;

  case INSIDE_ARRAY:

    if ( xmlStrcasecmp(name, (xmlChar *) "array") == 0 )
    {
      state = INSIDE_DATA;
    }

    break;

  default:
    ;
  }
}
//-----------------------------------------------------------------------------
void PXMLMesh::open(std::string filename)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
bool PXMLMesh::close()
{
  return state == DONE;
}
//-----------------------------------------------------------------------------
void PXMLMesh::readMesh(const xmlChar *name, const xmlChar **attrs)
{
  // Parse values
  std::string type = parseString(name, attrs, "celltype");
  uint gdim = parseUnsignedInt(name, attrs, "dim");
  
  // Create cell type to get topological dimension
  CellType* cell_type = CellType::create(type);
  uint tdim = cell_type->dim();
  delete cell_type;

  // Open mesh for editing
  editor.open(_mesh, CellType::string2type(type), tdim, gdim);
}
//-----------------------------------------------------------------------------
void PXMLMesh::readVertices(const xmlChar *name, const xmlChar **attrs)
{

  numParsedVertices = 0;
  numParsedCells = 0;

  // Parse values
  uint num_vertices = parseUnsignedInt(name, attrs, "size");

  uint pe_size = MPI::numProcesses();
  uint rank = MPI::processNumber();

  // Calculate a linear data distribution
  uint L = floor( (real) num_vertices / (real) pe_size);
  uint R = num_vertices % pe_size;
  numLocalVertices = (num_vertices +  pe_size - rank -1) / pe_size;

  startIndex_vert = rank * L + std::min(rank,R);
  endIndex_vert = startIndex_vert + ( numLocalVertices - 1);
  _mesh.distdata().set_global_numVertices(num_vertices);

  // Set number of vertices
  editor.initVertices(numLocalVertices);
}
//-----------------------------------------------------------------------------
void PXMLMesh::readCells(const xmlChar *name, const xmlChar **attrs)
{
  // Parse values
  //  uint num_cells = parseUnsignedInt(name, attrs, "size");

  numParsedCells = 0;

 //FIXME
  editor.initCells(1); 
  editor.close();
  MeshFunction<uint> pre_partition;
  _mesh.partition_geom(pre_partition); 
  _mesh.distribute(pre_partition);    

  for(VertexIterator vertex(_mesh); !vertex.end(); ++vertex)
    own_vertex[_mesh.distdata().get_global(*vertex)] = true;

}
//-----------------------------------------------------------------------------
void PXMLMesh::readVertex(const xmlChar *name, const xmlChar **attrs)
{
  // Read index
  uint v = parseUnsignedInt(name, attrs, "index");

  if(v < startIndex_vert || v > endIndex_vert)
    return;
  
  _mesh.distdata().set_map(numParsedVertices, v, 0);
  
  // Handle differently depending on geometric dimension
  switch ( _mesh.geometry().dim() )
    {
    case 1:
      {
	real x = parseReal(name, attrs, "x");
	editor.addVertex(numParsedVertices, x);
      }
      break;
    case 2:
      {
	real x = parseReal(name, attrs, "x");
	real y = parseReal(name, attrs, "y");
	editor.addVertex(numParsedVertices, x, y);
      }
      break;
    case 3:
      {
	real x = parseReal(name, attrs, "x");
	real y = parseReal(name, attrs, "y");
	real z = parseReal(name, attrs, "z");
	editor.addVertex(numParsedVertices, x, y, z);
      }
      break;
    default:
      error("Dimension of mesh must be 1, 2 or 3.");
    }
  numParsedVertices++;
}
//-----------------------------------------------------------------------------
void PXMLMesh::readInterval(const xmlChar *name, const xmlChar **attrs)
{
  // Check dimension
  if ( _mesh.topology().dim() != 1 )
    error("Mesh entity (interval) does not match dimension of mesh (%d).",
	  _mesh.topology().dim());
  
  // Parse values
  uint c  = parseUnsignedInt(name, attrs, "index");
  
  if(c < startIndex_cell || c > endIndex_cell)
    return;
  c = numParsedCells++;  
  
  uint v0 = parseUnsignedInt(name, attrs, "v0");
  uint v1 = parseUnsignedInt(name, attrs, "v1");
  // Add cell
  editor.addCell(c, v0, v1);
}
//-----------------------------------------------------------------------------
void PXMLMesh::readTriangle(const xmlChar *name, const xmlChar **attrs)
{
  // Check dimension
  if ( _mesh.topology().dim() != 2 )
    error("Mesh entity (triangle) does not match dimension of mesh (%d).",
		 _mesh.topology().dim());

  
// Parse values
//  uint c  = parseUnsignedInt(name, attrs, "index");

  uint v0 = parseUnsignedInt(name, attrs, "v0");
  uint v1 = parseUnsignedInt(name, attrs, "v1");
  uint v2 = parseUnsignedInt(name, attrs, "v2");

  // Return if no vertices are local
  if(!(own_vertex[v1] || own_vertex[v2] || own_vertex[v0]) || !own_vertex[v0])
    return;
  
  used_vertex[_mesh.distdata().get_local(v0, 0)] = true;
  if(own_vertex[v1])
    used_vertex[_mesh.distdata().get_local(v1, 0)] = true;
  if(own_vertex[v2])
    used_vertex[_mesh.distdata().get_local(v2, 0)] = true;

  // Add shared vertices to shared list
  if(!(own_vertex[v1] && own_vertex[v2] && own_vertex[v0])){
    if(!own_vertex[v1])
      shared_buffer[v1] = true;
    if(!own_vertex[v2])
      shared_buffer[v2] = true;
  }

  // Add cell to cell buffer
  cell_buffer.push_back(v0);
  cell_buffer.push_back(v1);
  cell_buffer.push_back(v2);
}
//-----------------------------------------------------------------------------
void PXMLMesh::readTetrahedron(const xmlChar *name, const xmlChar **attrs)
{
 // Check dimension
  if ( _mesh.topology().dim() != 3 )
    error("Mesh entity (tetrahedron) does not match dimension of mesh (%d).",
		 _mesh.topology().dim());

  // Parse values
  uint c  = parseUnsignedInt(name, attrs, "index");
  c = 0;
  uint v0 = parseUnsignedInt(name, attrs, "v0");
  uint v1 = parseUnsignedInt(name, attrs, "v1");
  uint v2 = parseUnsignedInt(name, attrs, "v2");
  uint v3 = parseUnsignedInt(name, attrs, "v3");
  
  // Return if no vertices are local
  if(!(own_vertex[v1] || own_vertex[v2] || own_vertex[v0] || own_vertex[v3])
     || !own_vertex[v0])
    return;

  used_vertex[_mesh.distdata().get_local(v0, 0)] = true;
  if(own_vertex[v1])
    used_vertex[_mesh.distdata().get_local(v1, 0)] = true;
  if(own_vertex[v2])
    used_vertex[_mesh.distdata().get_local(v2, 0)] = true;
  if(own_vertex[v3])
    used_vertex[_mesh.distdata().get_local(v3, 0)] = true;
  // Add shared vertices to shared list
  if(!(own_vertex[v1] && own_vertex[v2] && own_vertex[v0] && own_vertex[v3])){
    if(!own_vertex[v1])
      shared_buffer[v1] = true;
    if(!own_vertex[v2])
      shared_buffer[v2] = true;
    if(!own_vertex[v3])
      shared_buffer[v3] = true;
  }
  
  // Add cell to cell buffer
  cell_buffer.push_back(v0);
  cell_buffer.push_back(v1);
  cell_buffer.push_back(v2);
  cell_buffer.push_back(v3);  
  
}
//-----------------------------------------------------------------------------
void PXMLMesh::readMeshFunction(const xmlChar* name, const xmlChar** attrs)
{
  // Parse values
  const std::string id = parseString(name, attrs, "name");
  const std::string type = parseString(name, attrs, "type");
  const uint dim = parseUnsignedInt(name, attrs, "dim");
  const uint size = parseUnsignedInt(name, attrs, "size");

  // Only uint supported at this point
  if (strcmp(type.c_str(), "uint") != 0)
    error("Only uint-valued mesh data is currently supported.");

  // Check size
  _mesh.init(dim);
  if (_mesh.size(dim) != size)
    error("Wrong number of values for MeshFunction, expecting %d.", _mesh.size(dim));

  // Register data
  f = _mesh.data().createMeshFunction(id);
  dolfin_assert(f);
  f->init(_mesh, dim);

  // Set all values to zero
  *f = 0;
}
void PXMLMesh::readArray(const xmlChar* name, const xmlChar** attrs)
{
  // Parse values
  const std::string id = parseString(name, attrs, "name");
  const std::string type = parseString(name, attrs, "type");
  const uint size = parseUnsignedInt(name, attrs, "size");

  // Only uint supported at this point
  if (strcmp(type.c_str(), "uint") != 0)
    error("Only uint-valued mesh data is currently supported.");

  // Register data
  a = _mesh.data().createArray(id, size);
  dolfin_assert(a);
}
//-----------------------------------------------------------------------------
void PXMLMesh::readMeshEntity(const xmlChar* name, const xmlChar** attrs)
{
  // Read index
  const uint index = parseUnsignedInt(name, attrs, "index");

  // Read and set value
  dolfin_assert(f);
  dolfin_assert(index < f->size());
  const uint value = parseUnsignedInt(name, attrs, "value");
  f->set(index, value);
}
//-----------------------------------------------------------------------------
void PXMLMesh::readArrayElement(const xmlChar* name, const xmlChar** attrs)
{
  // Read index
  const uint index = parseUnsignedInt(name, attrs, "index");

  // Read and set value
  dolfin_assert(a);
  dolfin_assert(index < a->size());
  const uint value = parseUnsignedInt(name, attrs, "value");
  (*a)[index] = value;
}
//-----------------------------------------------------------------------------
void PXMLMesh::closeMesh()
{
  Mesh new_mesh;
  editor.open(new_mesh, _mesh.type().cellType(),
	      _mesh.topology().dim(), _mesh.geometry().dim());

  int rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  Array<uint> send_buff, send_buff_indices, send_buff_orphan;
  // Make room for own vertices and shared ones
  uint shared = 0;
  uint orphan = 0;
  for(uint i=0; i<shared_buffer.size(); i++)
    if(shared_buffer[i]){
      shared++;
      send_buff.push_back(i);
    }
   
  _map<uint,bool> assigned_orphan, ghost_vertex;
  for(uint i=0;i<_mesh.numVertices();i++){
    if(!used_vertex[i]){
      orphan++;
    }
  }

  int num_recv, max_nsh;
  Array<double> send_buff_coords;
  MPI_Status status;
  uint gdim = _mesh.geometry().dim();
  uint num_shared = shared;
  uint num_coords = gdim * num_shared;
  uint num_orphan = num_shared;
  real *shared_coords = new real[num_coords];      
  real *shp = &shared_coords[0];
  uint *shared_indices = new uint[num_shared];
  uint *shpi = &shared_indices[0];
  uint *shared_orphans = new uint[num_shared];
  uint *shpo = &shared_orphans[0];
  MPI_Allreduce(&num_shared, &max_nsh, 1, MPI_INT,MPI_MAX, MPI::DOLFIN_COMM); 
  uint *shvert = new uint[max_nsh];
  uint src,dest;
  _map<uint, uint> owner_map;
  // Exchange ghost points
  for(uint j = 1; j < pe_size ; j++){

    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_UNSIGNED, dest, 1,
		 shvert, max_nsh, MPI_UNSIGNED, src, 1, 
		 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&num_recv);

    for(int k=0; k<num_recv; k++)
      if(own_vertex[shvert[k]]){
	Vertex vertex(_mesh,_mesh.distdata().get_local(shvert[k], 0));
	send_buff_coords.push_back(vertex.point().x());
	send_buff_coords.push_back(vertex.point().y());
	if(gdim >2)
	  send_buff_coords.push_back(vertex.point().z());
	
	send_buff_indices.push_back(shvert[k]);
	if(!used_vertex[vertex.index()] &&
	   !assigned_orphan[vertex.index()]){
	  owner_map[vertex.index()] =  status.MPI_SOURCE;
	  send_buff_orphan.push_back(pe_size);
	  assigned_orphan[vertex.index()] = true;
	}
	else {
	  if(owner_map.count(vertex.index()) == 0)
	    send_buff_orphan.push_back(rank);	 
	  else
	    send_buff_orphan.push_back(owner_map[vertex.index()]);	 
	}
	ghost_vertex[_mesh.distdata().get_local(shvert[k], 0)]  = true;
      }

    MPI_Sendrecv(&send_buff_indices[0], send_buff_indices.size(), MPI_UNSIGNED,
		 src, 1, shpi, num_shared, MPI_UNSIGNED, dest, 1, 
		 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&num_recv);	
    num_shared -=num_recv;
    shpi += num_recv;

    MPI_Sendrecv(&send_buff_coords[0], send_buff_coords.size(), MPI_DOUBLE,
		 src, 2, shp, num_coords, MPI_DOUBLE, dest, 2, MPI::DOLFIN_COMM,
		 &status);
    MPI_Get_count(&status,MPI_DOUBLE,&num_recv);
    num_coords -=num_recv;
    shp +=num_recv;	

    MPI_Sendrecv(&send_buff_orphan[0], send_buff_orphan.size(), MPI_UNSIGNED,
		 src, 3, shpo, num_orphan, MPI_UNSIGNED, dest, 3, 
		 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status,MPI_UNSIGNED,&num_recv);
    num_orphan -= num_recv;
    shpo += num_recv;    

    send_buff_coords.clear();
    send_buff_indices.clear();
    send_buff_orphan.clear();

  }
  
  // Init new mesh
  editor.initVertices(_mesh.numVertices() + shared - orphan);
  //  new_mesh.distdata().init(_mesh.numVertices() + shared - orphan);
  new_mesh.distdata().set_global_numVertices(_mesh.distdata().global_numVertices());

  uint v=0;
  for(VertexIterator vertex(_mesh); !vertex.end(); ++vertex){    
    if(used_vertex[vertex->index()]){
      editor.addVertex(v,vertex->point());
      new_mesh.distdata().set_map(v, _mesh.distdata().get_global(*vertex), 0);
      if(ghost_vertex[vertex->index()])
	new_mesh.distdata().set_shared(v, 0);
      v++;
    }
  }

  //Add shared ghost vertices
  uint ci=0;
  for(uint i=0;i<shared;i++){
    new_mesh.distdata().set_map(v, shared_indices[i], 0);
    if(shared_orphans[i] < pe_size) { // Why...ugly hack to set ghost owner
      new_mesh.distdata().set_ghost(v, 0);
      new_mesh.distdata().set_ghost_owner(v, shared_orphans[i], 0);
    }
    switch(gdim)
      {
      case 2:
	editor.addVertex(v++,shared_coords[ci],shared_coords[ci+1]); break;
      case 3:
	editor.addVertex(v++,shared_coords[ci],shared_coords[ci+1],
			 shared_coords[ci+2]); 
	break;
	
      }
    ci += gdim;
  }  

  uint ndims =  _mesh.type().numVertices(_mesh.topology().dim());
  editor.initCells(cell_buffer.size() / ndims);
  uint c=0;
  for(uint i=0; i<cell_buffer.size(); i +=ndims){
    switch(ndims)
      {
      case 2:
	editor.addCell(c++,new_mesh.distdata().get_local(cell_buffer[i], 0),
		       new_mesh.distdata().get_local(cell_buffer[i+1], 0));
	break;	
      case 3:
	editor.addCell(c++,new_mesh.distdata().get_local(cell_buffer[i], 0),
		       new_mesh.distdata().get_local(cell_buffer[i+1], 0),
		       new_mesh.distdata().get_local(cell_buffer[i+2], 0));
	break;	
      case 4:
	editor.addCell(c++,new_mesh.distdata().get_local(cell_buffer[i], 0),
		       new_mesh.distdata().get_local(cell_buffer[i+1], 0), 
		       new_mesh.distdata().get_local(cell_buffer[i+2], 0),
		       new_mesh.distdata().get_local(cell_buffer[i+3], 0));
	break;	
      }
  }

  editor.close();
  _mesh = new_mesh;

  ghost_vertex.clear();
  assigned_orphan.clear();
  shared_buffer.clear();  
  send_buff.clear();
  cell_buffer.clear();
  used_vertex.clear();
  own_vertex.clear();
  delete[] shared_indices;
  delete[] shared_coords;
  delete[] shared_orphans;
  delete[] shvert;

}
//-----------------------------------------------------------------------------
#else
void PXMLMesh::startElement(const xmlChar* name, const xmlChar** attrs)
{
}
//-----------------------------------------------------------------------------
void PXMLMesh::endElement(const xmlChar* name)
{
}
//-----------------------------------------------------------------------------
void PXMLMesh::open(std::string filename)
{
}
//-----------------------------------------------------------------------------
bool PXMLMesh::close()
{
  return false;
}
//-----------------------------------------------------------------------------
#endif
