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
void MeshSmoother::maph0(Mesh& mesh, Mesh& sub, MeshFunction<int>& cell_map,
                         MeshFunction<real>& h0, MeshFunction<real>& subh0)
{
  subh0.init(sub, sub.topology().dim());

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;
    if (cell_map.get(cell) != -1)
    {
      subh0.set(cell_map(cell), h0.get(cell));
    }
  }
}

//-----------------------------------------------------------------------------
bool MeshSmoother::onBoundary(Cell& cell)
{
  int d = cell.dim();

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
void MeshSmoother::worstElement(Mesh& mesh, int& index,
                                MeshFunction<bool>& masked_cells)
{
  int d = mesh.topology().dim();
  mesh.init(d - 1, d);

  MeshQuality mqual(mesh);

  real mu_min = 1.0e12;
  index = -1;

  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    real qual = mqual.mean_ratio(cell);
    if (qual < mu_min && !onBoundary(cell) && !masked_cells.get(cell))
    {
      index = cell.index();
      mu_min = qual;
    }
  }
}

//-----------------------------------------------------------------------------
void MeshSmoother::elementNhood(Mesh& mesh, Cell& element,
                                MeshFunction<bool>& elements, int depth)
{
  elements.set(element, true);

  if (depth == 0) return;

  for (CellIterator c(element); !c.end(); ++c)
  {
    Cell& cell = *c;

    //elements.set(cell, true);
    elementNhood(mesh, cell, elements, depth - 1);
  }
}

//-----------------------------------------------------------------------------
void MeshSmoother::submesh(Mesh& mesh, Mesh& sub,
                           MeshFunction<bool>& smoothed_cells,
                           MeshFunction<int>& old2new_vertex,
                           MeshFunction<int>& old2new_cell)
{

  //dolfin_debug("Entering create submesh");

  old2new_vertex.init(mesh, 0);
  old2new_cell.init(mesh, mesh.topology().dim());

  int ncells = 0;
  int nvertices = 0;

  // Count cells and vertices in submesh
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    if (smoothed_cells.get(cell) == true)
    {
      ncells++;
    }
  }

  for (VertexIterator n(mesh); !n.end(); ++n)
  {
    Vertex& vertex = *n;

    bool included = false;

    for (CellIterator c(vertex); !c.end(); ++c)
    {
      Cell& cell = *c;

      if (smoothed_cells.get(cell) == true)
      {
        included = true;
      }
    }

    if (included)
    {
      nvertices++;
    }
  }

  // Get cell type
  const CellType& cell_type = mesh.type();

  unsigned int current_vertex = 0;
  unsigned int current_cell = 0;


  MeshDistributedData distdata(sub.topology().dim());
  MeshEditor editor(sub, cell_type.cellType(), mesh.geometry().dim());

  // Specify number of vertices and cells
  editor.init_vertices(nvertices);
  editor.init_cells(ncells);

  for (VertexIterator n(mesh); !n.end(); ++n)
  {
    Vertex& vertex = *n;

    bool included = false;

    for (CellIterator c(vertex); !c.end(); ++c)
    {
      Cell& cell = *c;

      if (smoothed_cells.get(cell) == true)
      {
        included = true;
      }
    }

    if (included)
    {
      old2new_vertex.set(vertex.index(), current_vertex);
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
      old2new_vertex.set(vertex.index(), -1);
    }
  }

  Array<unsigned int> cell_vertices(cell_type.num_entities(0));
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    if (smoothed_cells.get(cell) == true)
    {
      int cv_idx = 0;
      for (VertexIterator n(cell); !n.end(); ++n)
      {
        int id = old2new_vertex.get(n->index());
        if (id == -1)
        {
          cout << "broken: " << n->index() << endl;
        }
        cell_vertices[cv_idx++] = id;
      }

      old2new_cell.set(cell.index(), current_cell);
      distdata[c->dim()].set_map(current_cell, c->global_index());
      editor.add_cell(current_cell++, &cell_vertices[0]);

    }
    else
    {
      old2new_cell.set(cell.index(), -1);
    }
  }

  editor.close();
  sub.distdata() = distdata;
//  sub.distdata().set_invalid_numbering();
  sub.renumber();
}

}
