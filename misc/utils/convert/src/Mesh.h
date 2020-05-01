/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <njansson@kth.se> wrote this file. As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return Niclas Jansson
 * ----------------------------------------------------------------------------
 */

#ifndef DOLFIN_CONVERT_MESH_H
#define DOLFIN_CONVERT_MESH_H

struct Mesh
{

  Mesh() :
    gdim(0),
    cell_type(0),
    num_cells(0),
    num_vertices(0),
    cells(NULL),
    vertices(NULL) {}
  
  ~Mesh()
  {
    gdim = cell_type = num_cells = num_vertices = 0;
    
    if (cells)
      delete[] cells;
    
    if(vertices)
      delete[] vertices;
  }

  virtual void load_mesh(std::string& filename) = 0;

  uint32_t gdim;
  uint32_t cell_type;
  
  uint32_t num_cells;
  uint32_t num_vertices;
  
  uint32_t *cells;
  double *vertices;

};

#endif

