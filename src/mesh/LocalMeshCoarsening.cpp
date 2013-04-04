// Copyright (C) 2006 Johan Hoffman.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Balthasar Reuter, 2013.
// 
// First added:  2006-11-01
// Last changed: 2013-04-03

#include <dolfin/mesh/LocalMeshCoarsening.h>
#include <dolfin/mesh/CoarseningManager.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/TriangleCell.h>

#include <limits>

using namespace dolfin;
//-----------------------------------------------------------------------------
static inline real distance(real const * x0, real const * x1, uint dim)
{
  real sqrlength(0);

  if ( dim == 2 )
  {
    sqrlength = ( x0[0] - x1[0] )*( x0[0] - x1[0] ) + 
                ( x0[1] - x1[1] )*( x0[1] - x1[1] );
  }
  else if ( dim == 3 )
  {
    sqrlength = ( x0[0] - x1[0] )*( x0[0] - x1[0] ) + 
                ( x0[1] - x1[1] )*( x0[1] - x1[1] ) + 
                ( x0[2] - x1[2] )*( x0[2] - x1[2] );
  }
  else
    error("Unknown geometrical dimension!");

  return sqrt(sqrlength);
}
//-----------------------------------------------------------------------------
void LocalMeshCoarsening::coarsenMeshByEdgeCollapse(Mesh& mesh,
                                                    MeshFunction<bool>& cell_marker,
                                                    bool coarsen_boundary )
{
  uint init_num_cells = mesh.numCells();
  uint init_num_verts = mesh.numVertices();

  // check size of cell_marker
  if ( cell_marker.size() != mesh.numCells() )
    error( "Wrong dimension of cell_marker" );

  // Generate cell - edge connectivity if not yet generated
  mesh.init(mesh.topology().dim(), 1);

  // Generate edge - vertex connectivity if not yet generated
  mesh.init(1,0);

  // Create new mesh
  Mesh coarse_mesh(mesh);

  // Instantiate coarsening manager
  CoarseningManager manager(cell_marker, coarsen_boundary);

  uint num_cells_to_coarsen( manager.cells_to_coarsen().size() );
  message("%d cells selected for coarsening", num_cells_to_coarsen);

  // Coarsen until nothing happens anymore
  std::pair<bool,bool> result;
  uint prev_num_cells;
  do {
    prev_num_cells = mesh.numCells();

    // Try all the cells that are marked for coarsening
    for ( List<uint>::iterator c_it(manager.cells_to_coarsen().begin()) ; 
          c_it != manager.cells_to_coarsen().end() ; /* do nothing */ )
    {
      // try to coarsen the cell
      result = coarsenCell(mesh, coarse_mesh, manager, *c_it);

      // Coarsening successful: first result is true
      if ( result.first )
      {
        // remove cell from coarsening list
        c_it = manager.cells_to_coarsen().erase(c_it);

        // mesh changed: second result is true -> commit changes
        if ( result.second )
        {
          manager.updateDistdata(mesh, coarse_mesh);
          mesh = coarse_mesh; // TODO: swap instead of assignment possible??

          manager.vertex_map().setNewFineFromCoarseSize(mesh.numVertices());
          manager.vertex_map().commit();

          manager.cell_map().setNewFineFromCoarseSize(mesh.numCells());
          manager.cell_map().commit();
        }
      }
      else
      {
        // not successful: undo changes
        manager.vertex_map().revert();
        manager.cell_map().revert();

        // try the next one
        ++c_it;
      }
    }
  } while ( manager.migrate(mesh, prev_num_cells > mesh.numCells()) );

  message(
    "Mesh coarsening done. Deleted %d vertices and %d cells, %d vertices and %d cells remain",
     init_num_verts - mesh.numVertices(), init_num_cells - mesh.numCells(), 
     mesh.numVertices(), mesh.numCells() );
}
//-----------------------------------------------------------------------------
int LocalMeshCoarsening::selectEdge(Cell& c, CoarseningManager& manager)
{
  real lmin(std::numeric_limits<real>::max());
  int shortest_edge_index(-1);
  for ( EdgeIterator e_it(c) ; !e_it.end() ; ++e_it )
  {
    uint * verts = e_it->entities(0);
    real l = e_it->length();

    if ( 
        lmin < l &&                           // no shorter edge found before
        !(          // edge cannot be coarsened if both vertices are forbidden
          manager.isForbiddenVertex(
            manager.vertex_map().getFineFromCoarse(verts[0]) ) && 
          manager.isForbiddenVertex(
            manager.vertex_map().getFineFromCoarse(verts[1]) ) 
    ) ) 
    {
      lmin = l;
      shortest_edge_index = e_it->index();
    }
  }

  return shortest_edge_index;
}
//-----------------------------------------------------------------------------
bool LocalMeshCoarsening::selectEdge(Cell& c, CoarseningManager& manager, 
                                     uint * vertices)
{
  real lmin(std::numeric_limits<real>::max());
  bool edge_found(false);
  for ( VertexIterator v_it1(c) ; !v_it1.end() ; ++v_it1 )
  {
    VertexIterator v_it2(c);
    while( v_it2.pos() != v_it1.pos() ) ++v_it2;

    for ( ++v_it2 ; !v_it2.end() ; ++v_it2 )
    {
      uint dim = c.mesh().geometry().dim();
      real l = distance(v_it1->x(), v_it2->x(), dim);
      if ( 
          lmin > l &&                            // no shorter edge found before
        !(            // edge cannot be coarsened if both vertices are forbidden
          manager.isForbiddenVertex(
            manager.vertex_map().getFineFromCoarse(v_it1->index()) ) && 
          manager.isForbiddenVertex(
            manager.vertex_map().getFineFromCoarse(v_it2->index()) ) 
      ) ) 
      {
        lmin = l;
        vertices[0] = v_it1->index();
        vertices[1] = v_it2->index();
        edge_found = true;
      }
    }
  }
  return edge_found;
}
//-----------------------------------------------------------------------------
bool LocalMeshCoarsening::selectVertex(Edge& e, CoarseningManager& manager,
                                       uint& vertD, uint& vertR)
{
  uint * verts = e.entities(0);

  // Both end vertices forbidden: collapse not possible. Should not happen
  // since edge should not have been selected in the first place
  dolfin_assert( !(
    manager.isForbiddenVertex(
      manager.vertex_map().getFineFromCoarse(verts[0]) ) &&
    manager.isForbiddenVertex(
      manager.vertex_map().getFineFromCoarse(verts[1]) ) 
  ) );

  // verts[0] allowed and verts[1] forbidden: select verts[0] for collapse
  if ( !manager.isForbiddenVertex(
          manager.vertex_map().getFineFromCoarse(verts[0]) ) &&
        manager.isForbiddenVertex(
          manager.vertex_map().getFineFromCoarse(verts[1]) ) )
  {
    vertD = verts[0];
    vertR = verts[1];
  }
  // verts[0] forbidden and verts[1] allowed: select verts[1] for collapse
  else if ( manager.isForbiddenVertex(
              manager.vertex_map().getFineFromCoarse(verts[0]) ) &&
           !manager.isForbiddenVertex(
              manager.vertex_map().getFineFromCoarse(verts[1]) ) )
  {
    vertD = verts[1];
    vertR = verts[0];
  }
  // both allowed
  else
  {
    /*// verts[0] on a boundary and verts[1] not on a boundary: select verts[1]
    if (
       manager.isBoundaryVertex(
        manager.vertex_map().getFineFromCoarse(verts[0]) ) &&
      !manager.isBoundaryVertex(
        manager.vertex_map().getFineFromCoarse(verts[1]) ) )
    {
      vertD = verts[1];
      vertR = verts[0];
    }
    // verts[0] not on a boundary and verts[1] on a boundary: select verts[0]
    else if (
      !manager.isBoundaryVertex(
        manager.vertex_map().getFineFromCoarse(verts[0]) ) &&
       manager.isBoundaryVertex(
        manager.vertex_map().getFineFromCoarse(verts[1]) ) )
    {
      vertD = verts[0];
      vertR = verts[1];
    }
    // none on a boundary
    else*/
    {
      if ( verts[0] < verts[1] )
      {
        vertD = verts[1];
        vertR = verts[0];
      }
      else
      {
        vertD = verts[0];
        vertR = verts[1];
      }
    }
  }

  // Check if selected vertex is on a process boundary. In this case neighboring
  // entities have to be requested first from other processes
  if ( manager.isInteriorBoundaryVertex( 
          manager.vertex_map().getFineFromCoarse(vertD) ) )
  {
    manager.vertices_to_request().push_back( 
          manager.vertex_map().getFineFromCoarse(vertD) );
    return false;
  }
  else
  {
    return true;
  }
}
//-----------------------------------------------------------------------------
int LocalMeshCoarsening::selectVertex(uint * vertices, 
                                      CoarseningManager& manager)
{
  // Both end vertices forbidden: collapse not possible. Should not happen
  // since edge should not have been selected in the first place
  dolfin_assert( !(
    manager.isForbiddenVertex(
      manager.vertex_map().getFineFromCoarse(vertices[0]) ) &&
    manager.isForbiddenVertex(
      manager.vertex_map().getFineFromCoarse(vertices[1]) ) 
  ) );

  int ret(-1);

  // 0 allowed and 1 forbidden: select 0 for collapse
  if ( !manager.isForbiddenVertex(
          manager.vertex_map().getFineFromCoarse(vertices[0]) ) &&
        manager.isForbiddenVertex(
          manager.vertex_map().getFineFromCoarse(vertices[1]) ) )
  {
    ret = 0;
  }
  // 0 forbidden and 1 allowed: select 1 for collapse
  else if ( manager.isForbiddenVertex(
              manager.vertex_map().getFineFromCoarse(vertices[0]) ) &&
           !manager.isForbiddenVertex(
              manager.vertex_map().getFineFromCoarse(vertices[1]) ) )
  {
    ret = 1;
  }
  // both allowed: choose higher index for collapse
  else
  {
    if ( vertices[0] < vertices[1] )
      ret = 1;
    else
      ret = 0;
  }

  // Check if selected vertex is on a process boundary. In this case neighboring
  // entities have to be requested first from other processes
  if ( manager.isInteriorBoundaryVertex( 
          manager.vertex_map().getFineFromCoarse(vertices[ret]) ) )
  {
    manager.vertices_to_request().push_back( 
          manager.vertex_map().getFineFromCoarse(vertices[ret]) );
    return -1;
  }
  else
  {
    return ret;
  }
}
//-----------------------------------------------------------------------------
void LocalMeshCoarsening::regenerateCells(Mesh const & mesh, 
                                          MeshEditor& editor, 
                                          Vertex& vertex_to_remove, 
                                          uint vertR, uint c_id, 
                                          MeshFunction<bool> const & cells_to_remove, 
                                          CoarseningManager& manager)
{
  Array<uint> cell_vertices(mesh.type().numEntities(0));

  // iterate over all cells adjacent to removed vertex
  for ( CellIterator c_it(vertex_to_remove) ; !c_it.end() ; ++c_it )
  {
    // skip cells that have been removed
    if ( cells_to_remove.get(*c_it) ) continue;

    // collect vertices for new cell
    uint cv_id(0);
    for ( VertexIterator v_it(*c_it) ; !v_it.end() ; ++v_it, ++cv_id )
    {
      if ( v_it->index() == vertex_to_remove.index() )
        cell_vertices[cv_id] = manager.vertex_map().getNewCoarseFromCoarse(vertR);
      else
        cell_vertices[cv_id] = manager.vertex_map().getNewCoarseFromCoarse(v_it->index());
    }

    // add new cell
    editor.addCell(c_id, cell_vertices);

    // Update cell map
    manager.cell_map().setNew(manager.cell_map().getFineFromCoarse(c_it->index()), c_id);
    ++c_id;
  }
}
//-----------------------------------------------------------------------------
bool LocalMeshCoarsening::checkMesh(Vertex& removed_vertex, Mesh& coarse_mesh,
                                    CoarseningManager& manager)
{
  real vol_tol = 1.e-5;//1.e-3;

  // Check for inverted cells and new cell volumes of cells adjacent to 
  // removed vertex
  for ( CellIterator c_it(removed_vertex) ; !c_it.end() ; ++c_it ) 
  {
    int c_id( manager.cell_map().getNewCoarseFromCoarse(c_it->index()) );

    // consider only existing cells
    if ( c_id >= 0 )
    {
      Cell c(coarse_mesh, c_id);

      // check qm of new cell
      real qm = c.volume() / c.diameter();
      if ( qm < vol_tol )
      {
        //warning("Cell quality too low, qm = %f", qm);
        return false;
      }

      // check orientation of new cell
      if ( c_it->orientation() != c.orientation() )
      {
        //warning("Cell orientation inverted");
        return false;
      }
    }
  }

  return true;
}
//-----------------------------------------------------------------------------
std::pair<bool,bool> LocalMeshCoarsening::coarsenCell(Mesh& mesh, Mesh& coarse_mesh, 
                                                       CoarseningManager& manager,
                                                       uint cell_to_coarsen_id)
{
  // Check if cell has already been deleted
  if ( manager.cell_map().getCoarseFromFine(cell_to_coarsen_id) < 0 )
    return std::make_pair(true,false);
  Cell cell_to_coarsen(mesh, manager.cell_map().getCoarseFromFine(cell_to_coarsen_id));

//  // Select edge for collapse
//  int e_id = selectEdge(cell_to_coarsen, manager);
//
//  // if cell cannot be coarsened: simply return
//  if ( e_id < 0 ) return std::make_pair(true,false);
  uint verts[2];
  if ( !selectEdge(cell_to_coarsen, manager, verts) )
    return std::make_pair(true,false);

//  // select vertex for collapse
//  Edge edge_to_collapse(mesh,e_id);
//  uint vertD, vertR;
//  if ( !selectVertex(edge_to_collapse, manager, vertD, vertR) )
  int vert_idx = selectVertex( verts, manager );
  if ( vert_idx < 0 )
  {
    // Cannot be coarsened due to missing entities from other processes
    manager.cells_to_request().push_back( cell_to_coarsen_id );
    return std::make_pair(true,false);
  }
  uint vertD = verts[vert_idx];
  uint vertR = verts[!vert_idx];

  Vertex vertex_to_remove(mesh,vertD);
  Vertex vertex_to_keep(mesh,vertR);

  // Cells to remove: cells adjacent to collapsed edge
  MeshFunction<bool> cells_to_remove(mesh, mesh.topology().dim());
  cells_to_remove = false;
  uint num_cells_to_remove(0);
//  for ( CellIterator c_it(edge_to_collapse) ; !c_it.end() ; ++c_it )
//  {
//    cells_to_remove.set(c_it->index(),true);
//    ++num_cells_to_remove;

//    // Update cell map: cells don't exist anymore
//    manager.cell_map().setNew(
//      manager.cell_map().getFineFromCoarse(c_it->index()), -1 );
//  }
  for ( CellIterator c_it(vertex_to_remove) ; !c_it.end() ; ++c_it )
  {
    if ( c_it->incident(vertex_to_keep) )
    {
      cells_to_remove.set(c_it->index(),true);
      ++num_cells_to_remove;

      // Update cell map: cells don't exist anymore
      manager.cell_map().setNew(
      manager.cell_map().getFineFromCoarse(c_it->index()), -1 );
    }
  }

  // Cells to regenerate: cells adjacent to removed vertex excluding removed cells
  MeshFunction<bool> cells_to_regenerate(mesh, mesh.topology().dim());
  cells_to_regenerate = false;
  for ( CellIterator c_it(vertex_to_remove) ; !c_it.end() ; ++c_it )
  {
    if ( !cells_to_remove.get(c_it->index()) )
    {
      cells_to_regenerate.set(c_it->index(), true);

      // Update cell map: will be set to new values later
      manager.cell_map().setNew(
        manager.cell_map().getFineFromCoarse(c_it->index()), -1 );
    }
  }

  // MeshEditor for new mesh
  MeshEditor editor;
  editor.open(coarse_mesh, mesh.type().cellType(), 
              mesh.topology().dim(), mesh.geometry().dim());
  editor.initVertices(mesh.numVertices() - 1);
  editor.initCells(mesh.numCells() - num_cells_to_remove);

  // Add old vertices
  uint v_id(0);
  for ( VertexIterator v_it(mesh) ; !v_it.end() ; ++v_it )
  {
    if ( vertD == v_it->index() )
    {
      // update vertex map: vertex doesn't exist anymore
      manager.vertex_map().setNew(
        manager.vertex_map().getFineFromCoarse(vertD), -1 );
    }
    else
    {
      // add vertex to new mesh
      editor.addVertex(v_id, v_it->point());

      // update vertex map: vertex is now at position v_id
      manager.vertex_map().setNew(
        manager.vertex_map().getFineFromCoarse(v_it->index()), v_id );
      ++v_id;
    }
  }

  // Add old cells
  Array<uint> cell_vertices(mesh.type().numEntities(0));
  uint c_id(0);
  for ( CellIterator c_it(mesh) ; !c_it.end() ; ++c_it )
  {
    if ( !cells_to_remove.get(c_it->index()) && 
         !cells_to_regenerate.get(c_it->index()) )
    {
      // Build list of vertices
      uint cv_id(0);
      for ( VertexIterator v_it(*c_it) ; !v_it.end() ; ++v_it )
        cell_vertices[cv_id++] = manager.vertex_map().getNewCoarseFromCoarse(v_it->index());

      // add cell to new mesh
      editor.addCell(c_id, cell_vertices);

      // update cell map: cell is now at position c_id
      manager.cell_map().setNew(manager.cell_map().getFineFromCoarse(c_it->index()), c_id);
      ++c_id;
    }
  }

  // Add new cells
  regenerateCells(mesh, editor, vertex_to_remove, vertR, c_id, 
                  cells_to_remove, manager);

  // Finish editing
  editor.close();

  // Check quality
  return std::make_pair(
          checkMesh(vertex_to_remove, coarse_mesh, manager),
          true
         );
}
//-----------------------------------------------------------------------------