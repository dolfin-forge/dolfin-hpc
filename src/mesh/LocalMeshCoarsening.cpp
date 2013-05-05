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
#ifdef ____USE_D_MESH____
#include <dolfin/mesh/DMesh.h>
#include <dolfin/mesh/DCell.h>
#include <dolfin/mesh/DVertex.h>
#endif

#include <limits>
#include <algorithm>

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
#ifdef ____USE_D_MESH____
//-----------------------------------------------------------------------------
void LocalMeshCoarsening::coarsenMeshByEdgeCollapse(Mesh& mesh,
                                                    MeshFunction<bool>& cell_marker,
                                                    bool coarsen_boundary )
{
  uint init_num_cells = mesh.numCells();
  uint init_num_verts = mesh.numVertices();

  dolfin_assert( &(cell_marker.mesh()) == &mesh );

  // check size of cell_marker
  if ( cell_marker.size() != mesh.numCells() )
    error( "Wrong dimension of cell_marker" );

  // Instantiate coarsening manager
  CoarseningManager manager;
  manager.init(cell_marker, coarsen_boundary);

  dolfin_assert( manager.checkDCellNumbering(mesh.numCells() - 1) );

  uint num_cells_to_coarsen( manager.cells_to_coarsen().size() );
  message("[%d] Mesh coarsening", MPI::processNumber());

  // Coarsen until nothing happens anymore
  uint prev_num_cells_coarsened, num_cells_coarsened(0);
  int result;
  do {
    prev_num_cells_coarsened = num_cells_coarsened;

    // Try all the cells that are marked for coarsening
    for ( List<DCell*>::iterator c_it(manager.cells_to_coarsen().begin()) ; 
          c_it != manager.cells_to_coarsen().end() ; /* do nothing */ )
    {
      // try to coarsen the cell
      result = coarsenCell(manager, *c_it);

      // Coarsening not successful: try the next one
      if ( result < 0 ) 
        ++c_it;
      else // successful: remove cell from coarsening list
      {
        c_it = manager.cells_to_coarsen().erase(c_it);
        num_cells_coarsened += result;
      }
    }
  } while ( manager.migrate(mesh, prev_num_cells_coarsened < num_cells_coarsened) );

  Mesh omesh;
  manager.dmesh()->exp(omesh);
  mesh = omesh;

  message("[%d] Mesh coarsening done.", MPI::processNumber());
}
//-----------------------------------------------------------------------------
bool LocalMeshCoarsening::selectEdge(DCell* c, CoarseningManager& manager, 
                                     DVertex * vertices[])
{
  real lmin(std::numeric_limits<real>::max());
  bool edge_found(false);
  for ( std::vector<DVertex *>::iterator v_it1(c->vertices.begin()) ; 
        v_it1 != c->vertices.end() ; ++v_it1 )
  {
    DVertex * v1 = *v_it1;
    std::vector<DVertex *>::iterator v_it2 = v_it1;

    for ( ++v_it2 ; v_it2 != c->vertices.end() ; ++v_it2 )
    {
      DVertex * v2 = *v_it2;
      real l = v1->p.distance(v2->p);
      if ( 
          lmin > l &&                            // no shorter edge found before
        !(            // edge cannot be coarsened if both vertices are forbidden
          manager.isForbiddenVertex( v1->id ) && 
          manager.isForbiddenVertex( v2->id )
      ) ) 
      {
        lmin = l;
        vertices[0] = v1;
        vertices[1] = v2;
        edge_found = true;
      }
    }
  }
  return edge_found;
}
//-----------------------------------------------------------------------------
int LocalMeshCoarsening::selectVertex(DVertex * vertices[], 
                                      CoarseningManager& manager)
{
  // Both end vertices forbidden: collapse not possible. Should not happen
  // since edge should not have been selected in the first place
  dolfin_assert( !( manager.isForbiddenVertex( vertices[0]->id ) &&
                    manager.isForbiddenVertex( vertices[1]->id ) ) );

  int ret(-1);

  // 0 allowed and 1 forbidden: select 0 for collapse
  if ( !manager.isForbiddenVertex( vertices[0]->id ) &&
        manager.isForbiddenVertex( vertices[1]->id ) )
  {
    ret = 0;
  }
  // 0 forbidden and 1 allowed: select 1 for collapse
  else if ( manager.isForbiddenVertex( vertices[0]->id ) &&
           !manager.isForbiddenVertex( vertices[1]->id ) )
  {
    ret = 1;
  }
  // both allowed: choose higher index for collapse
  else
  {
    if ( vertices[0]->id < vertices[1]->id )
      ret = 1;
    else
      ret = 0;
  }

  // Check if selected vertex is on a process boundary. In this case neighboring
  // entities have to be requested first from other processes
  if ( manager.isInteriorBoundaryVertex( vertices[ret]->id ) )
  {
    manager.vertices_to_request().push_back( vertices[ret]->id );
    return -1;
  }
  else
  {
    return ret;
  }
}
//-----------------------------------------------------------------------------
bool LocalMeshCoarsening::checkMesh(std::list<DCell *>& cells_to_regenerate,
                                    std::vector<uint>& cells_to_regenerate_orient)
{
  real vol_tol = 1.e-5;//1.e-3;

  // Check for inverted cells and new cell volumes of cells adjacent to 
  // removed vertex
  std::vector<uint>::iterator o_it(cells_to_regenerate_orient.begin());
  for ( std::list<DCell *>::iterator c_it(cells_to_regenerate.begin()) ; 
        c_it != cells_to_regenerate.end() ; ++c_it, ++o_it ) 
  {
    DCell * dc = *c_it;
    dolfin_assert( !dc->deleted );

    // check qm of new cell
    real qm = dc->volume() / dc->diameter();
    if ( qm < vol_tol )
    {
      //warning("Cell quality too low, qm = %f", qm);
      return false;
    }

    // check orientation of new cell
    if ( dc->orientation() != *o_it )
    {
      //warning("Cell orientation inverted");
      return false;
    }
  }

  return true;
}
//-----------------------------------------------------------------------------
int LocalMeshCoarsening::coarsenCell(CoarseningManager& manager,
                                     DCell * cell_to_coarsen)
{
  // Check if cell has already been deleted
  if ( cell_to_coarsen->deleted )
    return 0;

  // Select edge for collapse
  DVertex * verts[2];
  if ( !selectEdge(cell_to_coarsen, manager, verts) )
    return 0;

  int vert_idx = selectVertex( verts, manager );
  if ( vert_idx < 0 )
  {
    // Cannot be coarsened due to missing entities from other processes
    manager.cells_to_request().push_back( cell_to_coarsen->id );
    return -1;
  }
  DVertex * vertex_to_remove = verts[vert_idx];
  DVertex * vertex_to_keep = verts[!vert_idx];

  // Cells to remove: all cells containing both vertices,
  // can be found as intersection of sorted cell lists of the two vertices
  std::list<DCell *> vertex_to_remove_cells(vertex_to_remove->cells);
  std::list<DCell *> vertex_to_keep_cells(vertex_to_keep->cells);
  vertex_to_remove_cells.sort();
  vertex_to_keep_cells.sort();

  std::list<DCell *> cells_to_remove;
  std::set_intersection( vertex_to_remove_cells.begin(), 
                         vertex_to_remove_cells.end(),
                         vertex_to_keep_cells.begin(), 
                         vertex_to_keep_cells.end(),
                         std::back_inserter(cells_to_remove) );

  // Cells to regenerate: all cells adjacent to removed vertex, that are not 
  // marked for removal
  // can be found as difference of sorted cell lists of vertex and removal
  std::list<DCell *> cells_to_regenerate;
  std::set_difference( vertex_to_remove_cells.begin(), 
                       vertex_to_remove_cells.end(),
                       cells_to_remove.begin(),
                       cells_to_remove.end(),
                       std::back_inserter(cells_to_regenerate) );

  dolfin_assert( cells_to_regenerate.size() + cells_to_remove.size() == 
                  vertex_to_remove_cells.size() );

  // Save cell orientations for checkMesh
  std::vector<uint> cells_to_regenerate_orientations(cells_to_regenerate.size());
  std::vector<uint>::iterator o_it(cells_to_regenerate_orientations.begin());
  for ( std::list<DCell *>::iterator c_it(cells_to_regenerate.begin()) ; 
        c_it != cells_to_regenerate.end() ; ++c_it, ++o_it )
    *o_it = (*c_it)->orientation();

  // Mark vertex as removed in global list
  manager.dmesh()->removeVertex(vertex_to_remove);

  // Delete cells from cell-list in vertices and mark as removed in global list
  for ( std::list<DCell *>::iterator c_it(cells_to_remove.begin()) ;
        c_it != cells_to_remove.end() ; ++c_it )
  {
    DCell * dc = *c_it;
    for ( std::vector<DVertex *>::iterator v_it(dc->vertices.begin()) ;
          v_it != dc->vertices.end() ; ++v_it )
    {
      DVertex * dv = *v_it;
      std::list<DCell *>::iterator rm_it = std::find( dv->cells.begin(),
                                                      dv->cells.end(), dc );
      dolfin_assert( rm_it != dv->cells.end() );
      dv->cells.erase(rm_it);
    }
    manager.dmesh()->removeCell(*c_it);
  }

  // Regenerate cells, i. e. replace vertex in vertices-vector of cell and add
  // to cell-vector of kept vertex
  for ( std::list<DCell *>::iterator c_it(cells_to_regenerate.begin()) ; 
        c_it != cells_to_regenerate.end() ; ++c_it )
  {
    DCell * dc = *c_it;
    std::vector<DVertex *>::iterator v_it = std::find(dc->vertices.begin(),
                                                     dc->vertices.end(),
                                                     vertex_to_remove);
    *v_it = vertex_to_keep;
    vertex_to_keep->cells.push_back(dc);
  }

  // Check quality
  if ( checkMesh(cells_to_regenerate, cells_to_regenerate_orientations) )
    return cells_to_remove.size();

  // Quality not ok: revert changes (i. e. unmark vertex and cells as deleted,
  // add cells to local lists in vertices and undo regeneration)
  vertex_to_remove->deleted = false;
  
  for ( std::list<DCell *>::iterator c_it(cells_to_remove.begin()) ;
        c_it != cells_to_remove.end() ; ++c_it )
  {
    DCell * dc = *c_it;
    for ( std::vector<DVertex *>::iterator v_it(dc->vertices.begin()) ;
          v_it != dc->vertices.end() ; ++v_it )
    {
      DVertex * dv = *v_it;
      dv->cells.push_back(dc);
    }
    dc->deleted = false;
  }

  for ( std::list<DCell *>::iterator c_it(cells_to_regenerate.begin()) ; 
        c_it != cells_to_regenerate.end() ; ++c_it )
  {
    DCell * dc = *c_it;
    
    std::list<DCell *>::iterator rm_it = std::find(vertex_to_keep->cells.begin(),
                                                   vertex_to_keep->cells.end(),
                                                   dc);
    dolfin_assert( rm_it != vertex_to_keep->cells.end() );
    vertex_to_keep->cells.erase(rm_it);

    std::vector<DVertex *>::iterator v_it = std::find(dc->vertices.begin(),
                                                     dc->vertices.end(),
                                                     vertex_to_keep);
    *v_it = vertex_to_remove;
  }

  return -1;
}
//-----------------------------------------------------------------------------
#else // ____USE_D_MESH____
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

  // Create new mesh
  Mesh coarse_mesh;//(mesh);

  // Instantiate coarsening manager
  CoarseningManager manager;
  manager.init(cell_marker, coarsen_boundary);

  uint num_cells_to_coarsen( manager.cells_to_coarsen().size() );
  message("[%d] Mesh coarsening", MPI::processNumber());

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

  if (MPI::numProcesses() > 1) 
  {
    //mesh.distdata().invalid_numbering();
    //mesh.distdata().invalid_ownership();
    mesh.renumber();
  }

  message("[%d] Mesh coarsening done.", MPI::processNumber() );
}
//-----------------------------------------------------------------------------
#ifdef ____AVOID_TOPOLOGY_INIT____
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
#else // ____AVOID_TOPOLOGY_INIT____
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
        lmin > l &&                           // no shorter edge found before
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
#endif // ____AVOID_TOPOLOGY_INIT____
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

#ifdef ____AVOID_TOPOLOGY_INIT____
  uint verts[2];
  if ( !selectEdge(cell_to_coarsen, manager, verts) )
    return std::make_pair(true,false);
#else
  // Select edge for collapse
  int e_id = selectEdge(cell_to_coarsen, manager);

  // if cell cannot be coarsened: simply return
  if ( e_id < 0 ) 
    return std::make_pair(true,false);
#endif

  uint vertD, vertR;

#ifdef ____AVOID_TOPOLOGY_INIT____
  int vert_idx = selectVertex( verts, manager );
  if ( vert_idx < 0 )
  {
    // Cannot be coarsened due to missing entities from other processes
    manager.cells_to_request().push_back( cell_to_coarsen_id );
    return std::make_pair(true,false);
  }
  vertD = verts[vert_idx];
  vertR = verts[!vert_idx];
#else
  // select vertex for collapse
  Edge edge_to_collapse(mesh,e_id);
  if ( !selectVertex(edge_to_collapse, manager, vertD, vertR) )
  {
    // Cannot be coarsened due to missing entities from other processes
    manager.cells_to_request().push_back( cell_to_coarsen_id );
    return std::make_pair(true,false);
  }
#endif

  Vertex vertex_to_remove(mesh,vertD);
  Vertex vertex_to_keep(mesh,vertR);

  // Cells to remove: cells adjacent to collapsed edge
  MeshFunction<bool> cells_to_remove(mesh, mesh.topology().dim());
  cells_to_remove = false;
  uint num_cells_to_remove(0);

#ifdef ____AVOID_TOPOLOGY_INIT____
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
#else
  for ( CellIterator c_it(edge_to_collapse) ; !c_it.end() ; ++c_it )
  {
    cells_to_remove.set(c_it->index(),true);
    ++num_cells_to_remove;

    // Update cell map: cells don't exist anymore
    manager.cell_map().setNew(
      manager.cell_map().getFineFromCoarse(c_it->index()), -1 );
  }
#endif

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
#endif // ____USE_D_MESH____
//-----------------------------------------------------------------------------