// Copyright (C) 2009-2011 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First  added: 2009
// Last changed: 2011-02-22

#include <fstream>
#include <dolfin/common/types.h>
#include <dolfin/la/Vector.h>
#include <dolfin/io/BinaryFile.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>

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

  message(1, "Saved vector  to file %s in binary format.", filename.c_str());  
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
  std::ofstream fp(filename.c_str(), std::ofstream::binary);

  // Write Header
  int dim = mesh.geometry().dim();
  fp.write((char *)&dim, sizeof(int));

  int type = 0;
  CellType::Type cell_type = mesh.type().cellType();
  if (CellType::type2string(cell_type).c_str() == "tetrahedron")
    type = 1;
  fp.write((char *)&type, sizeof(int));

  // Write vertices 
  for (VertexIterator v(mesh); !v.end(); ++v)
    fp.write((char *)v->x(), dim * sizeof(double));

  // Write cells
  for (CellIterator c(mesh); !c.end(); ++c)
    fp.write((char *)c->entities(0), 4 * sizeof(int));
  

  
}
//----------------------------------------------------------------------------

