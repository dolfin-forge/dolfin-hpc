// Copyright (C) 2009-2011 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First  added: 2009
// Last changed: 2011-04-02

#include <fstream>
#include <dolfin/common/types.h>
#include <dolfin/la/Vector.h>
#include <dolfin/io/BinaryFile.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>



#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace dolfin;

//----------------------------------------------------------------------------
BinaryFile::BinaryFile(const std::string filename) : GenericFile(filename)
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

  std::ifstream fp(filename.c_str(), std::ifstream::binary);

  uint size;
  fp.read((char *)&size, sizeof(uint));
  
  real *values = new real[size];  

  fp.read((char *)values, size * sizeof(real));
  fp.close();

  x.init(size);
  x.set(values);
  delete[] values;  
}
//----------------------------------------------------------------------------
void BinaryFile::operator<<(GenericVector& x)
{
  
  std::ofstream fp(filename.c_str(), std::ofstream::binary);

  real *values = new real[x.local_size()];
  uint size = x.local_size();
  x.get(values);

  fp.write((char *)&size, sizeof(uint));
  fp.write((char *)values, x.local_size() * sizeof(real));
  fp.close();

  delete[] values;     

  message(1, "Saved vector to file %s in binary format.", filename.c_str());  
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
    MPI_File fh;
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
		  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);    


    int dim, type, num_vertices;
    MPI_File_read_all(fh, &dim, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_all(fh, &type, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_all(fh, &num_vertices, 1, MPI_INT, MPI_STATUS_IGNORE);
    message("dim: %d type: %d num_vertice: %d", dim, type, num_vertices);
    uint pe_size = MPI::numProcesses();
    uint pe_rank = MPI::processNumber();
    
    uint L = floor( (real) num_vertices / (real) pe_size);
    uint R = num_vertices % pe_size;
    uint local_vertices = (num_vertices + pe_size - pe_rank -1 ) / pe_size;

    
    uint offset = 0;
    uint vertex_data = dim * local_vertices;
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&vertex_data, &offset, 1, 
	       MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
    MPI_Scan(&vertex_data, &offset, 1, 
	     MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
    offset -= vertex_data;
#endif

    real *vertex_buffer = new real[vertex_data];
    MPI_File_read_at(fh, 3*sizeof(int) + offset * sizeof(double),
		     vertex_buffer, vertex_data, MPI_DOUBLE, MPI_STATUS_IGNORE);
              
    int num_cells;
    MPI_File_read_at(fh, 3*sizeof(int) + dim * num_vertices * sizeof(double),
		     &num_cells, 1, MPI_INT, MPI_STATUS_IGNORE);
    message("Numcells: %d", num_cells);
    uint local_cells = (num_cells + pe_size - pe_rank - 1 ) / pe_size;    

    offset = 0;
    uint cell_data = (3 + type) * local_cells;
#if ( MPI_VERSION > 1 )
    MPI_Exscan(&cell_data, &offset, 1, 
	       MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
    MPI_Scan(&cell_data, &offset, 1, 
	     MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
    offset -= cell_data;
#endif

    int *cell_buffer = new int[cell_data];
    MPI_File_read_at(fh, 4 * sizeof(int) + dim * num_vertices * sizeof(double) 
		     + offset * sizeof(int), cell_buffer, cell_data, 
		     MPI_INT, MPI_STATUS_IGNORE);

    message("Pre parse cells");
    // Parse cells
    std::vector<atomic_cell> cells;
    std::vector<uint> *non_local_cells = new std::vector<uint>[pe_size];
    std::vector<uint> *ghosts = new std::vector<uint>[pe_size];
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
       	cells.push_back(cell);
	
	if(vertex_owner(L, R, cell.v2) != pe_rank)
	  ghosts[vertex_owner(L, R, cell.v2)].push_back(cell.v2);

	if(vertex_owner(L, R, cell.v3) != pe_rank)
	  ghosts[vertex_owner(L, R, cell.v3)].push_back(cell.v3);

	if (type == 1)
	  if(vertex_owner(L, R, cell.v4) != pe_rank)
	    ghosts[vertex_owner(L, R, cell.v4)].push_back(cell.v4);
      }
      else
      {
	non_local_cells[vertex_owner(L, R, cell_buffer[i])].push_back(cell.v1);
	non_local_cells[vertex_owner(L, R, cell_buffer[i])].push_back(cell.v2);
	non_local_cells[vertex_owner(L, R, cell_buffer[i])].push_back(cell.v3);
	non_local_cells[vertex_owner(L, R, cell_buffer[i])].push_back(cell.v4);
      }

    }


    message("Pre communication");
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
	
	cells.push_back(cell);

	if(vertex_owner(L, R, cell.v2) != pe_rank)
	  ghosts[vertex_owner(L, R, cell.v2)].push_back(cell.v2);
	
	if(vertex_owner(L, R, cell.v3) != pe_rank)
	  ghosts[vertex_owner(L, R, cell.v3)].push_back(cell.v3);
	
	if (type == 1)
	  if(vertex_owner(L, R, cell.v4) != pe_rank)
	    ghosts[vertex_owner(L, R, cell.v4)].push_back(cell.v4);       	
      }      
    }

    delete[] cell_buffer;
    


    delete[] recv_buffer;
    delete[] vertex_buffer;

    
    MPI_File_close(&fh);
    message("MPI I/O: Done reading file");
  }
}  
//----------------------------------------------------------------------------
void BinaryFile::operator<<(Mesh& mesh)
{

  int dim = mesh.geometry().dim();

  int type = 0;
  CellType::Type cell_type = mesh.type().cellType();
  if (CellType::type2string(cell_type).c_str() == "tetrahedron")
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

    MPI_File fh;

    /* FIXME:
     * Add MPI_Info data 
     * Split and cleanup implementation
     */
    
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
		  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    // Write Header
    MPI_File_write_all(fh, &dim, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_write_all(fh, &type, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_write_all(fh, &num_vertices, 1, MPI_INT, MPI_STATUS_IGNORE);

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

    MPI_File_write_at(fh, 3*sizeof(int) + offset*sizeof(double), 
		      vertex_buffer, vertex_data, 
		      MPI_DOUBLE, MPI_STATUS_IGNORE);

    delete[] vertex_buffer;

    // Write Cells
    MPI_File_write_at_all(fh, 3*sizeof(int) + 
			  (dim * num_vertices) * sizeof(double),
			  &num_cells, 1, MPI_INT, MPI_STATUS_IGNORE);


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

    MPI_File_write_at(fh, 4*sizeof(int) + 
		      (dim * num_vertices) * sizeof(double) + 
		      offset * sizeof(int), 
		      cell_buffer,local_cell_entities , 
		      MPI_INT, MPI_STATUS_IGNORE);    

    delete[] cell_buffer;


    MPI_File_close(&fh);
    
  }

  message(1, "Saved mesh to file %s in binary format.", filename.c_str());    
}
//----------------------------------------------------------------------------


