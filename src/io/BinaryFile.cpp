// Copyright (C) 2009-2011 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First  added: 2009
// Last changed: 2011-04-01

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
    mesh.renumber();

    MPI_File fh;

    /* FIXME:
     * Add MPI_Info data 
     */
    
    MPI_File_open(dolfin::MPI::DOLFIN_COMM, (char *) filename.c_str(),
		  MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    // Write Header
    if (MPI::processNumber() == 0) 
    {
      MPI_File_write(fh, &dim, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write(fh, &type, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_write(fh, &num_vertices, 1, MPI_INT, MPI_STATUS_IGNORE);
    }

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
    if (MPI::processNumber() == 0)
      MPI_File_write_at(fh, 3*sizeof(int) + 
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

