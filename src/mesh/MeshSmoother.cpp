// Copyright (C) 2014 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// Refactoring of classes from UNICORN.
//
// First added:  2014-06-12
// Last changed: 2014-06-12

#include <dolfin/mesh/MeshSmoother.h>

#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshQuality.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshSmoother::MeshSmoother(Mesh& mesh) :
    MeshDependent(mesh)
{
}

//-----------------------------------------------------------------------------
MeshSmoother::~MeshSmoother()
{
}

//-----------------------------------------------------------------------------
void MeshSmoother::maph0(MeshValues<int, Cell>& cell_map,
                         MeshValues<real, Cell>& h0,
                         MeshValues<real, Cell>& subh0)
{
  for (CellIterator c(cell_map.mesh()); !c.end(); ++c)
  {
    if (cell_map(*c) != -1) { subh0(cell_map(*c)) = h0(*c); }
  }
}

//-----------------------------------------------------------------------------
bool MeshSmoother::onBoundary(Cell& cell)
{
  uint const d = cell.dim();
  for (FacetIterator f(cell); !f.end(); ++f)
  {
    if (f->num_entities(d) == 1)
    {
      return true;
    }
  }
  return false;
}

//-----------------------------------------------------------------------------
void MeshSmoother::worstElement(int& index, MeshValues<bool, Cell>& masked_cells)
{
  Mesh& mesh = masked_cells.mesh();
  MeshQuality mqual(mesh);

  real mu_min = 1.0e12;
  index = -1;

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    real qual = mqual.mean_ratio(*c);
    if (qual < mu_min && !onBoundary(*c) && !masked_cells(*c))
    {
      index = c->index();
      mu_min = qual;
    }
  }
}

//-----------------------------------------------------------------------------
void MeshSmoother::elementNhood(Cell& element,
                                MeshValues<bool, Cell>& elements, int depth)
{
  elements(element) = true;

  if (depth == 0) return;

  for (CellIterator c(element); !c.end(); ++c)
  {
    elementNhood(*c, elements, depth - 1);
  }
}

//-----------------------------------------------------------------------------
void MeshSmoother::submesh(Mesh& sub,
                           MeshValues<bool, Cell>& smoothed_cells,
                           MeshValues<int, Vertex>& old2new_vertex,
                           MeshValues<int, Cell>& old2new_cell)
{
  int ncells = 0;
  int nvertices = 0;

  // Count cells and vertices in submesh
  Mesh& mesh = smoothed_cells.mesh();
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    if (smoothed_cells(*c) == true) { ++ncells; }
  }

  bool included = false;
  for (VertexIterator n(smoothed_cells.mesh()); !n.end(); ++n,
       nvertices+= included, included = false)
  {
    for (CellIterator c(*n); !c.end(); ++c)
    {
      if (smoothed_cells(*c) == true) { included = true; }
    }
  }

  // Get cell type
  const CellType& cell_type = mesh.type();

  unsigned int current_vertex = 0;
  unsigned int current_cell = 0;

  MeshEditor editor(sub, cell_type.cellType(), mesh.geometry_dimension());

  // Specify number of vertices and cells
  editor.init_vertices(nvertices);
  editor.init_cells(ncells);

  MeshDistributedData& distdata = sub.distdata();

  for (VertexIterator n(mesh); !n.end(); ++n)
  {
    Vertex& vertex = *n;

    bool included = false;

    for (CellIterator c(vertex); !c.end(); ++c)
    {
      Cell& cell = *c;

      if (smoothed_cells(cell) == true)
      {
        included = true;
      }
    }

    if (included)
    {
      old2new_vertex(vertex.index()) = current_vertex;
      editor.add_vertex(current_vertex, vertex.x());
      distdata[0].set_map(current_vertex,
                       mesh.distdata()[0].get_global(vertex.index()));
      if (mesh.distdata()[0].is_ghost(vertex.index()))
      {
        distdata[0].set_ghost(current_vertex, vertex.owner());
      }

      current_vertex++;
    }
    else
    {
      old2new_vertex(vertex.index()) = -1;
    }
  }

  Array<uint> cell_vertices(cell_type.num_entities(0));
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    if (smoothed_cells(cell) == true)
    {
      int cv_idx = 0;
      for (VertexIterator n(cell); !n.end(); ++n)
      {
        int id = old2new_vertex(n->index());
        if (id == -1)
        {
          cout << "broken: " << n->index() << endl;
        }
        cell_vertices[cv_idx++] = id;
      }

      old2new_cell(cell.index()) = current_cell;
      distdata[c->dim()].set_map(current_cell, c->global_index());
      editor.add_cell(current_cell++, &cell_vertices[0]);

    }
    else
    {
      old2new_cell(cell.index()) = -1;
    }
  }

  editor.close();
}

}
