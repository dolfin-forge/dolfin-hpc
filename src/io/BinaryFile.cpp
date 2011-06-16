// Copyright (C) 2009-2011 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First  added: 2009
// Last changed: 2011-06-16

#include <algorithm>
#include <cstring>
#include <fstream>
#include <dolfin/common/types.h>
#include <dolfin/la/Vector.h>
#include <dolfin/io/BinaryFile.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/dolfin.h>


#ifdef ENABLE_MPIIO
#include <mpi.h>
#endif

using namespace dolfin;

//----------------------------------------------------------------------------
BinaryFile::BinaryFile(const std::string filename) : GenericFile(filename), 
						     _t(0)
{
  type = "Binary";
}
//----------------------------------------------------------------------------
BinaryFile::BinaryFile(const std::string filename, real &t) : 
  GenericFile(filename), _t(&t)
{
  type = "Binary";
}
//----------------------------------------------------------------------------
BinaryFile::~BinaryFile()
{
  // Do nothing
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(GenericVector& x)
{

  uint size;

#ifdef ENABLE_MPIIO
  uint offset[2];
  uint pe_rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  
  MPI_File fh;
  MPI_Offset byte_offset;
  BinaryFileHeader hdr;
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
		MPI_MODE_RDONLY , MPI_INFO_NULL, &fh);
  MPI_File_read_all(fh, &hdr, sizeof(BinaryFileHeader), 
		    MPI_BYTE, MPI_STATUS_IGNORE);
  
  hdr_check(hdr, BINARY_VECTOR_DATA, pe_size);

  byte_offset = sizeof(BinaryFileHeader);
  MPI_File_read_at_all(fh, byte_offset + pe_rank * 2 * sizeof(uint),
		       &offset[0], 2, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  size = offset[1];
  byte_offset += pe_size * 2 * sizeof(uint);

#else
  std::ifstream fp(filename.c_str(), std::ifstream::binary);
  fp.read((char *)&size, sizeof(uint));
#endif 

  real *values = new real[size];  

#ifdef ENABLE_MPIIO
  MPI_File_read_at_all(fh, byte_offset + offset[0] * sizeof(real),
		   values, offset[1], MPI_DOUBLE, MPI_STATUS_IGNORE);
  MPI_File_close(&fh);
#else
  fp.read((char *)values, size * sizeof(real));
  fp.close();
#endif

  x.init(size);
  x.set(values);
  delete[] values;  
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(GenericVector& x)
{
  uint size = x.local_size();
  real *values = new real[size];

#ifdef ENABLE_MPIIO
  uint offset[2] = {0,0};
  offset[0] = x.offset();
  offset[1] = size;
#endif

  x.get(values);

#ifdef ENABLE_MPIIO  
  BinaryFileHeader hdr;
  uint pe_rank = MPI::processNumber();  
  hdr.magic = BINARY_MAGIC;
  hdr.pe_size = MPI::numProcesses();
  hdr.type = BINARY_VECTOR_DATA;  
#ifdef HAVE_BIG_ENDIAN
  hdr.bendian = 1;
#else
  hdr.bendian = 0;
#endif

  MPI_File fh;
  MPI_Offset byte_offset;
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
		MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

  MPI_File_write_all(fh, &hdr, sizeof(BinaryFileHeader) , 
		     MPI_BYTE, MPI_STATUS_IGNORE);
  byte_offset = sizeof(BinaryFileHeader);
  MPI_File_write_at_all(fh, byte_offset + pe_rank * 2 *sizeof(uint), 
		    &offset[0], 2, MPI_UNSIGNED, MPI_STATUS_IGNORE);
  byte_offset += hdr.pe_size * 2 * sizeof(uint);
  MPI_File_write_at_all(fh, byte_offset + offset[0]*sizeof(real),
		    values, offset[1], MPI_DOUBLE, MPI_STATUS_IGNORE);

  MPI_File_close(&fh);
#else
  std::ofstream fp(filename.c_str(), std::ofstream::binary);
  fp.write((char *)&size, sizeof(uint));
  fp.write((char *)values, x.local_size() * sizeof(real));
  fp.close();
#endif

  delete[] values;     
  
  message(1, "Saved vector to file %s in binary format.", filename.c_str());  
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(Function &u)
{
  std::pair<Function*, std::string> f(&u, "U");
  std::vector<std::pair<Function*, std::string> > tmp;
  tmp.push_back(f);
  write_function(tmp);
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(std::vector<std::pair<Function*, std::string> >& f) 
{
  write_function(f);
}
//----------------------------------------------------------------------------
void BinaryFile::write_function(std::vector<std::pair<Function*, 
				std::string> >& f) 
{

#ifdef ENABLE_MPIIO

  nameUpdate(counter);

  BinaryFileHeader hdr;
  uint pe_rank = MPI::processNumber();  
  hdr.magic = BINARY_MAGIC;
  hdr.pe_size = MPI::numProcesses();
  hdr.type = BINARY_FUNCTION_DATA;  
#ifdef HAVE_BIG_ENDIAN
  hdr.bendian = 1;
#else
  hdr.bendian = 0;
#endif

  MPI_File fh;
  MPI_Offset byte_offset;
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) bin_filename.c_str(),
		MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  MPI_File_write_all(fh, &hdr, sizeof(BinaryFileHeader) , 
		     MPI_BYTE, MPI_STATUS_IGNORE);
  byte_offset = sizeof(BinaryFileHeader);

  uint n_func = f.size();
  MPI_File_write_all(fh, &n_func, 1, MPI_INT, MPI_STATUS_IGNORE);
  byte_offset += sizeof(int);

  // Assume same mesh for all data arrays
  Mesh& mesh = f[0].first->mesh();

  for (std::vector<std::pair<Function*, std::string> >::iterator it = f.begin();
       it != f.end(); it++) 
  {
    Function* u = it->first;
    
    std::string& name = it->second;
    const uint rank = u->rank();
    if(rank > 1)
      error("Only scalar and vectors functions can be saved in Binary.");
    
    // Get number of components
    const uint dim = u->dim(0);
    
    // Allocate memory for function values at vertices
    uint size = mesh.numVertices();
    for (uint i = 0; i < u->rank(); i++)
      size *= u->dim(i);
    real *values = new real[size];
    u->vector().get(values);

    size = u->dim(0) * (mesh.numVertices() - mesh.distdata().num_ghost(0));      
    
    BinaryFunctionHeader f_hdr;
    f_hdr.dim = dim;
    f_hdr.size = u->dim(0) * mesh.distdata().global_numVertices();
    if (name.length() >FNAME_LENGTH)
      error("Function name too long.");      
    strcpy(&f_hdr.name[0], name.c_str());
    if (_t) 	
      f_hdr.t = *_t;
    else 
      f_hdr.t = counter;
    MPI_File_write_at_all(fh, byte_offset, &f_hdr, sizeof(BinaryFunctionHeader), 
			  MPI_BYTE, MPI_STATUS_IGNORE);
    byte_offset += sizeof(BinaryFunctionHeader);

    MPI_File_write_at_all(fh, byte_offset + u->vector().offset()*sizeof(real),
			  values, size, MPI_DOUBLE, MPI_STATUS_IGNORE);
    byte_offset += u->dim(0) * (mesh.distdata().global_numVertices() * sizeof(real));

    delete[] values;
  }
  
  MPI_File_close(&fh);
  
  counter++;
#else
  error("MPI I/O is required to save functions in Binary.");
#endif
}
//----------------------------------------------------------------------------
void BinaryFile::operator>>(Mesh& mesh)
{
 
  if (MPI::numProcesses() == 1) 
  {    
    std::ifstream fp(filename.c_str(), std::ifstream::binary);
    
    int celltype, gdim;
    fp.read((char *)&gdim, sizeof(int));  
    fp.read((char *)&celltype, sizeof(int));  

    std::string type;
    if (celltype == 0)
      type = "triangle";
    else if(celltype == 1)
      type = "tetrahedron";
    else
      error("Unknown Cell type");

    
    // Create cell type to get topological dimension
    CellType* cell_type = CellType::create(type);
    uint tdim = cell_type->dim();
    delete cell_type;
    
    // Open mesh for editing
    MeshEditor editor;
    editor.open(mesh, CellType::string2type(type), tdim, gdim);
    
    // Read vertex data
    int n_vertices;
    fp.read((char *)&n_vertices, sizeof(int));  
    editor.initVertices(n_vertices);
    
    real *vertex_data = new real[n_vertices * gdim];
    fp.read((char *)vertex_data, n_vertices * gdim * sizeof(double));
    
    int v = 0; 
    for (int i = 0; i < n_vertices * gdim; v++, i += gdim) 
    {
      switch(gdim)
      {
      case 2:
	editor.addVertex(v, vertex_data[i], vertex_data[i+1]); 
	break;
      case 3:
	editor.addVertex(v, vertex_data[i], 
			 vertex_data[i+1], vertex_data[i+2]);
	break;
      default:
	error("Dimension of mesh must be 1, 2 or 3.");      
      }
    }  
    delete[] vertex_data;
    
    // Read cell data
    int n_cells;
    fp.read((char *)&n_cells, sizeof(int));  
    editor.initCells(n_cells);
    int *cell_data = new int[n_cells * (3 + celltype)];
    fp.read((char *)cell_data, n_cells * (3 + celltype) * sizeof(int));
    
    int c = 0;
    for (int i = 0; i < n_cells * (3 + celltype); c++, i+= (3 + celltype)) 
    {
      switch(celltype)
      {
      case 0:
	editor.addCell(c, cell_data[i], cell_data[i+1], cell_data[i+2]);
	break;
      case 1:
	editor.addCell(c, cell_data[i], cell_data[i+1], 
		       cell_data[i+2], cell_data[i+3]);
	break;
      }
    }
    delete[] cell_data;
    editor.close();
    fp.close();
    
  }
  else
  {
#ifdef ENABLE_MPIIO
    
    uint pe_size = MPI::numProcesses();
    uint pe_rank = MPI::processNumber();

    MPI_File fh;
    MPI_Offset byte_offset;
    BinaryFileHeader hdr;
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
		  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);    


    int dim, type, num_vertices;
    MPI_File_read_all(fh, &hdr, sizeof(BinaryFileHeader), 
		      MPI_BYTE, MPI_STATUS_IGNORE);
    hdr_check(hdr, BINARY_MESH_DATA, pe_size);
    MPI_File_read_all(fh, &dim, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_all(fh, &type, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_all(fh, &num_vertices, 1, MPI_INT, MPI_STATUS_IGNORE);

    byte_offset = sizeof(BinaryFileHeader) + 3 * sizeof(int);

    uint L = floor( (real) num_vertices / (real) pe_size);
    uint R = num_vertices % pe_size;
    uint local_vertices = (num_vertices + pe_size - pe_rank -1 ) / pe_size;

    
    uint offset[2] = {0,0};
    uint vertex_data[2] = {local_vertices, dim * local_vertices};
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&vertex_data[0], &offset[0], 2, 
	       MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
    MPI_Scan(&vertex_data[0], &offset[0], 2, 
	     MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
    offset[0] -= vertex_data[0];
    offset[1] -= vertex_data[1];
#endif

    real *vertex_buffer = new real[vertex_data[1]];
    MPI_File_read_at_all(fh, byte_offset + offset[1] * sizeof(double), 
			 vertex_buffer, vertex_data[1], 
			 MPI_DOUBLE, MPI_STATUS_IGNORE);
    byte_offset += dim * num_vertices * sizeof(double);
              
    int num_cells;
    MPI_File_read_at_all(fh, byte_offset, &num_cells, 
			 1, MPI_INT, MPI_STATUS_IGNORE);
    byte_offset += sizeof(int);

    uint local_cells = (num_cells + pe_size - pe_rank - 1 ) / pe_size;    

    offset[1] = 0;
    uint cell_data = (3 + type) * local_cells;
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&cell_data, &offset[1], 1, 
	       MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
    MPI_Scan(&cell_data, &offset[1], 1, 
	     MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
    offset[1] -= cell_data;
#endif

    int *cell_buffer = new int[cell_data];
    MPI_File_read_at_all(fh, byte_offset + offset[1] * sizeof(int), 
			 cell_buffer, cell_data, MPI_INT, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);
    
    _set<uint> ghosted_entities;

    // Parse cells
    std::vector<atomic_cell> cells;
    std::vector<uint> *non_local_cells = new std::vector<uint>[pe_size];
    std::set<uint> used_vertices;

    atomic_cell cell;
    for (int i = 0; i < cell_data; i+= (3 + type)) 
    {

      cell.v1 = cell_buffer[i];
      cell.v2 = cell_buffer[i+1];
      cell.v3 = cell_buffer[i+2];
      if (type == 1)       
	cell.v4 = cell_buffer[i+3];      

      if (vertex_owner(L, R, cell_buffer[i]) == pe_rank)
      {
	used_vertices.insert(cell.v1);
       	cells.push_back(cell);
	
	if(vertex_owner(L, R, cell.v2) != pe_rank)
	  ghosted_entities.insert(cell.v2);
	else
	  used_vertices.insert(cell.v2);

	if(vertex_owner(L, R, cell.v3) != pe_rank)
	  ghosted_entities.insert(cell.v3);
	else
	  used_vertices.insert(cell.v3);

	if (type == 1)
	  if(vertex_owner(L, R, cell.v4) != pe_rank) 
	    ghosted_entities.insert(cell.v4);
	  else
	    used_vertices.insert(cell.v4);
      }
      else
      {
	non_local_cells[vertex_owner(L, R, cell_buffer[i])].push_back(cell.v1);
	non_local_cells[vertex_owner(L, R, cell_buffer[i])].push_back(cell.v2);
	non_local_cells[vertex_owner(L, R, cell_buffer[i])].push_back(cell.v3);
	if (type == 1)
	  non_local_cells[vertex_owner(L, R, cell_buffer[i])].push_back(cell.v4);
      }

    }

    /*
     * FIXME
     * Reduce communication in this section
     */
    uint local_max = 0;
    for (int i = 0; i < pe_size; i++) 
      local_max = std::max(local_max, (uint) non_local_cells[i].size());
      
    uint buff_size = 0;
    MPI_Allreduce(&local_max, &buff_size, 1, MPI_UNSIGNED, MPI_MAX, MPI::DOLFIN_COMM);


    uint *recv_buffer = new uint[buff_size];    

    // Exchange data
    MPI_Status status;
    int num_recv, src, dest;
    for (int i = 1; i < pe_size; i++) 
    {
      src = (pe_rank - i + pe_size) % pe_size;
      dest = (pe_rank + i) % pe_size;
      
      MPI_Sendrecv(&non_local_cells[dest][0], non_local_cells[dest].size(), 
		   MPI_UNSIGNED, dest, 1, recv_buffer, buff_size, 
		   MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &num_recv);
       
      // Add received cells
      for (int j = 0; j < num_recv; j += (3 + type)) 
      {
	cell.v1 = recv_buffer[j];
	cell.v2 = recv_buffer[j+1];
	cell.v3 = recv_buffer[j+2];
	if (type == 1)
	  cell.v4 = recv_buffer[j+3];
	
	used_vertices.insert(cell.v1);
	cells.push_back(cell);

	if(vertex_owner(L, R, cell.v2) != pe_rank) 
	  ghosted_entities.insert(cell.v2);
	else
	  used_vertices.insert(cell.v2);
	
	if(vertex_owner(L, R, cell.v3) != pe_rank) 
	  ghosted_entities.insert(cell.v3);
	else
	  used_vertices.insert(cell.v3);
	
	if (type == 1)
	  if(vertex_owner(L, R, cell.v4) != pe_rank) 
	    ghosted_entities.insert(cell.v4);
	  else
	    used_vertices.insert(cell.v4);
      }      
    }

    // Number of vertices in mesh, used + ghosts
    uint num_local_vertices = used_vertices.size() + ghosted_entities.size();    


    std::string celltype;
    if (type == 0)
      celltype = "triangle";
    else if(type == 1)
      celltype = "tetrahedron";
    else
      error("Unknown Cell type");
    
    // Create cell type to get topological dimension
    CellType* cell_type = CellType::create(celltype);
    uint tdim = cell_type->dim();
    delete cell_type;
    


    std::set<uint> all_vertices, orphaned_vertices;
    for(uint i = 0; i < local_vertices; i++) 
      all_vertices.insert(offset[0] + i);
    std::set_difference(all_vertices.begin(), all_vertices.end(),
			used_vertices.begin(), used_vertices.end(),
			std::inserter(orphaned_vertices, orphaned_vertices.end()));

    // Open mesh for editing
    MeshEditor editor;
    editor.open(mesh, CellType::string2type(celltype), tdim, dim);
    num_local_vertices = all_vertices.size() + ghosted_entities.size();
    editor.initVertices(num_local_vertices);
    editor.initCells(cells.size());
    

    MeshDistributedData distdata;

    /* Add local indices
     * Use start_index for global number
     *
     * Exchange ghost vertices
     * Use ghosts to know global number
     *
     * Add cells, remember to map global->local
     */
   
    uint local_vertex_index = 0;
    uint v_index; 
    for(std::set<uint>::iterator it = used_vertices.begin();
	it != used_vertices.end(); local_vertex_index++, it++) 
    {
      v_index = *it - offset[0];

      distdata.set_map(local_vertex_index, *it, 0);
            
      switch(dim)
      {
      case 2:
	editor.addVertex(local_vertex_index,
			 vertex_buffer[(dim * v_index)],
			 vertex_buffer[(dim * v_index) + 1]);
	break;
      case 3:
	editor.addVertex(local_vertex_index,
			 vertex_buffer[(dim * v_index)],
			 vertex_buffer[(dim * v_index) + 1],
			 vertex_buffer[(dim * v_index) + 2]);
	break;
      }
    }
    
    _map<uint, uint> new_owner;
    // Exchange ghost points
    std::vector<uint> *ghosts = new std::vector<uint>[pe_size];    
    for(_set<uint>::iterator it = ghosted_entities.begin(); 
	it != ghosted_entities.end(); it++) 
      ghosts[vertex_owner(L, R, *it)].push_back(*it);       	
    
    local_max = 0;
    for (int i = 0; i < pe_size; i++) 
      local_max = std::max(local_max, (uint) ghosts[i].size());

    buff_size = 0;
    MPI_Allreduce(&local_max, &buff_size, 1, MPI_UNSIGNED, MPI_MAX, MPI::DOLFIN_COMM);
 
    delete[] recv_buffer;
    recv_buffer = new uint[buff_size]; 
    
    uint *send_new_owner = new uint[buff_size];
    uint *recv_new_owner = new uint[buff_size];

    real *send_buffer_coords = new real[buff_size * dim];
    real *recv_buffer_coords = new real[buff_size * dim];

    for (int i = 1; i < pe_size; i++) 
    {
      src = (pe_rank - i + pe_size) % pe_size;
      dest = (pe_rank + i) % pe_size;
      
      MPI_Sendrecv(&ghosts[dest][0], ghosts[dest].size(), 
		   MPI_UNSIGNED, dest, 1, recv_buffer, buff_size, 
		   MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &num_recv);

      /*
       * Check if orphaned
       * Send back global number and orphaned info 
       */

      real *sp = &send_buffer_coords[0];
      uint *np = &send_new_owner[0];
      for (int k = 0; k < num_recv; k++) 
      {
	
	v_index = recv_buffer[k] - offset[0];   
	
	if (orphaned_vertices.find(recv_buffer[k]) != orphaned_vertices.end()) {
	  if (new_owner.find(recv_buffer[k]) == new_owner.end())
	    new_owner[recv_buffer[k]] = src;
	  *(np++) = new_owner[recv_buffer[k]];
	}
	else {
	  *(np++) = pe_rank;
	}

	for (uint l = 0; l < dim; l++)
	  *(sp++)= vertex_buffer[(dim * v_index) + l];
      }
      MPI_Sendrecv(send_new_owner, num_recv, MPI_UNSIGNED, src, 1, 
		   recv_new_owner, buff_size, MPI_UNSIGNED, dest, 1,
		   MPI::DOLFIN_COMM, &status);

      MPI_Sendrecv(send_buffer_coords, (num_recv * dim), MPI_DOUBLE, src, 1, 
		   recv_buffer_coords, (buff_size * dim), MPI_DOUBLE, dest, 1,
		   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &num_recv);
      
      int g_i = 0;
      for (int k = 0; k < num_recv; local_vertex_index++, g_i++, k += dim) 
      {
	distdata.set_map(local_vertex_index, ghosts[dest][g_i], 0);

	if (recv_new_owner[g_i] != pe_rank) {
	  distdata.set_ghost(local_vertex_index, 0);
	  distdata.set_ghost_owner(local_vertex_index, recv_new_owner[g_i], 0);	  
	}

	switch(dim)
	{
	case 2:
	editor.addVertex(local_vertex_index,
			 recv_buffer_coords[k], recv_buffer_coords[k+1]);
	break;
	case 3:
	editor.addVertex(local_vertex_index,
			 recv_buffer_coords[k], recv_buffer_coords[k+1], recv_buffer_coords[k+2]);
	break;
	}
	
      }      
    }

    uint local_cell_index = 0;
    for(std::vector<atomic_cell>::iterator it =  cells.begin();
	it != cells.end(); local_cell_index++, it++) 
    {
       switch(type)
      {
      case 0:
	editor.addCell(local_cell_index, distdata.get_local(it->v1, 0), 
		       distdata.get_local(it->v2, 0), distdata.get_local(it->v3, 0));
	break;
      case 1:
	editor.addCell(local_cell_index, distdata.get_local(it->v1, 0), distdata.get_local(it->v2, 0),
		       distdata.get_local(it->v3, 0), distdata.get_local(it->v4, 0));
	break;
      }
    }


    editor.close();
    mesh.distdata() = distdata;


    delete[] cell_buffer;
    delete[] recv_buffer;
    delete[] vertex_buffer;
    delete[] recv_buffer_coords;
    delete[] send_buffer_coords;
    delete[] recv_new_owner;
    delete[] send_new_owner;
    
    for (uint i = 0; i < pe_size; i++) 
    {
      non_local_cells[i].clear();
      ghosts[i].clear();
    }
    delete[] non_local_cells;
    delete[] ghosts;

#else
    error("MPI and MPI I/O is required for reading binary meshes in parallel");
#endif

  }
}  
//----------------------------------------------------------------------------
void BinaryFile::operator<<(Mesh& mesh)
{

  int dim = mesh.geometry().dim();

  int type = 0;
  CellType::Type cell_type = mesh.type().cellType();
  if (CellType::type2string(cell_type) == "tetrahedron")
    type = 1;

  int num_vertices = (MPI::numProcesses() > 1 ? 
		      mesh.distdata().global_numVertices() : 
		      mesh.numVertices());
  int num_cells = (MPI::numProcesses() > 1 ? 
		      mesh.distdata().global_numCells() : mesh.numCells());
  if (MPI::numProcesses() == 1) 
  {    
    std::ofstream fp(filename.c_str(), std::ofstream::binary);
    
    // Write Header
    fp.write((char *)&dim, sizeof(int));
    
    fp.write((char *)&type, sizeof(int));
    
    // Write vertices 
    fp.write((char *)&num_vertices, sizeof(int));      
    for (VertexIterator v(mesh); !v.end(); ++v)
      fp.write((char *)v->x(), dim * sizeof(real));
    
    // Write cells
    fp.write((char *) &num_cells, sizeof(int));   
    for (CellIterator c(mesh); !c.end(); ++c)
      fp.write((char *)c->entities(0), (3 + type) * sizeof(int));
        
    fp.close();
  }
  else 
  {

#ifdef ENABLE_MPIIO
    BinaryFileHeader hdr;
    hdr.magic = BINARY_MAGIC;
    hdr.pe_size = MPI::numProcesses();
    hdr.type = BINARY_MESH_DATA;
#ifdef HAVE_BIGENDIAN
    hdr.bendian = 1;
#else
    hdr.bendian = 0;
#endif

    
    MPI_File fh;
    MPI_Offset byte_offset;
    /* FIXME:
     * Add MPI_Info data 
     * Split and cleanup implementation
     */
    
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
		  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    // Write Header
    MPI_File_write_all(fh, &hdr, sizeof(BinaryFileHeader), MPI_BYTE, 
		       MPI_STATUS_IGNORE);
    MPI_File_write_all(fh, &dim, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_write_all(fh, &type, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_write_all(fh, &num_vertices, 1, MPI_INT, MPI_STATUS_IGNORE);
    byte_offset = sizeof(BinaryFileHeader) + 3 * sizeof(int);

    // Write vertices
    uint offset = 0;
    uint vertex_data = dim * (mesh.numVertices() - mesh.distdata().num_ghost(0));
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&vertex_data, &offset, 1, 
	       MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
    MPI_Scan(&vertex_data, &offset, 1, 
	     MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
    offset -= vertex_data;
#endif

    real *vertex_buffer = new real[vertex_data];
    real *vp = &vertex_buffer[0];
    for (VertexIterator v(mesh); !v.end(); ++v) 
    {
      if(!mesh.distdata().is_ghost(v->index(), 0))
      {
	*(vp++) = v->x()[0];
	*(vp++) = v->x()[1];
	if (dim == 3)
	  *(vp++) = v->x()[2];
      }
    }

    MPI_File_write_at_all(fh, byte_offset + offset*sizeof(double), 
		      vertex_buffer, vertex_data, 
		      MPI_DOUBLE, MPI_STATUS_IGNORE);
    byte_offset += dim * num_vertices * sizeof(double);
    delete[] vertex_buffer;

    // Write Cells
    MPI_File_write_at_all(fh, byte_offset,&num_cells, 
			  1, MPI_INT, MPI_STATUS_IGNORE);
    byte_offset += sizeof(int);

    offset = 0;
    uint local_cell_entities = (3+type) * mesh.numCells();
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&local_cell_entities, &offset, 1, 
	       MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
    MPI_Scan(&local_cell_entities, &offset, 1, 
	     MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
    offset -= local_cell_entities;
#endif
    
    int *cell_buffer = new int[local_cell_entities];
    int *cp = &cell_buffer[0];
    for (CellIterator c(mesh); !c.end(); ++c) 
      for (uint i = 0; i < c->numEntities(0); i++) 
	*(cp++) = mesh.distdata().get_global(c->entities(0)[i], 0);

    MPI_File_write_at_all(fh, byte_offset + offset * sizeof(int),
		      cell_buffer,local_cell_entities , 
		      MPI_INT, MPI_STATUS_IGNORE);    

    delete[] cell_buffer;


    MPI_File_close(&fh);
    
#else
    error("MPI and MPI I/O is required for writing a mesh in parallel");
#endif
    
  }

  message(1, "Saved mesh to file %s in binary format.", filename.c_str());    
}
//----------------------------------------------------------------------------
void BinaryFile::nameUpdate(const int counter)
{
  std::string filestart, extension;
  std::ostringstream fileid, newfilename;
  
  fileid.fill('0');
  fileid.width(6);
  
  filestart.assign(filename, 0, filename.find("."));
  extension.assign(filename, filename.find("."), filename.size());
  
  fileid << counter;
  newfilename << filestart << fileid.str() << ".bin";
  
  bin_filename = newfilename.str();
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(MeshFunction<int>& meshfunction)
{
  write_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(MeshFunction<unsigned int>& meshfunction)
{
  write_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(MeshFunction<double>& meshfunction)
{
  write_meshfunction(meshfunction);
}
//----------------------------------------------------------------------------
template<class T>
void BinaryFile::write_meshfunction(T& meshfunction)
{

#ifdef ENABLE_MPIIO
  nameUpdate(counter);

  Mesh& mesh  = meshfunction.mesh();
  T *values = new T[meshfunction.size()];
  T *vp = &values[0];

  BinaryFileHeader hdr;
  uint pe_rank = MPI::processNumber();  
  hdr.magic = BINARY_MAGIC;
  hdr.pe_size = MPI::numProcesses();
  hdr.type = BINARY_MESH_FUNCTION_DATA;  
#ifdef HAVE_BIG_ENDIAN
  hdr.bendian = 1;
#else
  hdr.bendian = 0;
#endif

  MPI_File fh;
  MPI_Offset byte_offset;
  MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) bin_filename.c_str(),
		MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

  MPI_File_write_all(fh, &hdr, sizeof(BinaryFileHeader) , 
		     MPI_BYTE, MPI_STATUS_IGNORE);
  byte_offset = sizeof(BinaryFileHeader);

  uint local_size = 0;
  uint mfunc_type = 0;
  if ( meshfunction.dim() == mesh.topology().dim()) 
  {        
    MPI_File_write_at_all(fh, byte_offset, &mfunc_type, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += sizeof(uint);

    for (CellIterator c(mesh); !c.end(); ++c)
      *(vp++) = meshfunction.get(c->index());

    local_size = mesh.numCells();
  }
  else if ( meshfunction.dim() == 0) 
  {
    mfunc_type = 1;
    MPI_File_write_at_all(fh, byte_offset, &mfunc_type, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
    byte_offset += sizeof(uint);

    for (VertexIterator v(mesh); !v.end(); ++v)
      if (!mesh.distdata().is_ghost(v->index(), 0))
	*(vp++) = meshfunction.get(v->index());
    
    local_size = mesh.numVertices() - mesh.distdata().num_ghost(0);
  }
  else
    error("Binary output of mesh functions is implemenetd for cell/vertex-based functions only.");    

  uint offset = 0;
#if ( MPI_VERSION > 1 )
  MPI_Exscan(&local_size, &offset, 1, 
	     MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
  MPI_Scan(&local_size, &offset, 1, 
	   MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
  offset -= local_size;
#endif

  MPI_File_write_at_all(fh, byte_offset + offset * sizeof(T), values,
			local_size * sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);

  MPI_File_close(&fh);
  
  counter++;

#else
  error("MPI I/O required for writing mesh functions to binary files");
#endif
  
}
//----------------------------------------------------------------------------


