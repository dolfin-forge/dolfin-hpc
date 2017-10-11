// Copyright (C) 2006 Johan Hoffman.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Balthasar Reuter, 2013.
//
// First added:  2006-11-01
// Last changed: 2013-06-01


#include <list>

#include <dolfin/math/basic.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/LocalMeshCoarsening.h>
#include <dolfin/mesh/CoarseningManager.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/TriangleCell.h>
#include <dolfin/mesh/DMesh.h>
#include <dolfin/mesh/DCell.h>
#include <dolfin/mesh/DVertex.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/mesh/LoadBalancer.h>
#include <dolfin/common/timing.h>

#include <limits>
#include <algorithm>

using namespace dolfin;

//-----------------------------------------------------------------------------
void LocalMeshCoarsening::coarsenMeshByEdgeCollapse(MeshValues<bool, Cell>& cell_marker,
                                                    bool coarsen_boundary )
{
  begin("Coarsening simplicial mesh by edge collapse.");

  // Instantiate coarsening manager
  CoarseningManager manager(cell_marker, coarsen_boundary);

  uint num_cells_to_coarsen( manager.cells_to_coarsen().size() );
  message("%d cells selected on process %d",
	  num_cells_to_coarsen, MPI::rank());

  // Coarsen until nothing happens anymore
  uint prev_num_cells_coarsened, num_cells_coarsened(0);
  int result;
  do {
    prev_num_cells_coarsened = num_cells_coarsened;

    // Try all the cells that are marked for coarsening
    for ( List< std::pair<DCell*,uint> >::iterator
          c_it(manager.cells_to_coarsen().begin()) ;
          c_it != manager.cells_to_coarsen().end() ;  )
    {
      // try to coarsen the cell
      result = coarsenCell(manager, c_it->first, c_it->second);

      // Coarsening not successful: try the next one
      if ( result < 0 )
      {
        ++(c_it->second);
        // give up on this cell if failed several times
        // attempt count is one larger than actual number of attempts for
        // simplicity during migration phase
        if ( c_it->second > 10 )
          c_it = manager.cells_to_coarsen().erase(c_it);
        else
          ++c_it;
      }
      else // successful: remove cell from coarsening list
      {
        c_it = manager.cells_to_coarsen().erase(c_it);
        num_cells_coarsened += result;
      }
    }
  } while ( manager.migrate(num_cells_coarsened - prev_num_cells_coarsened) );

  Mesh omesh;
  manager.dmesh()->exp(omesh);
  cell_marker.mesh().swap(omesh);

  end();
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
      real l = v1->p.dist(v2->p);
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
                                      CoarseningManager& manager,
                                      uint attempts)
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
  // both allowed: alternate the points between attempts to overcome previous
  // failures due to mesh quality
  else
  {
    if ( attempts % 2 == 0 && vertices[0]->id < vertices[1]->id )
      ret = 1;
    else
      ret = 0;
  }

  // Check if selected vertex is on a process boundary. In this case neighboring
  // entities have to be requested first from other processes
  //if ( manager.isInteriorBoundaryVertex( vertices[ret]->id ) )
  if ( vertices[ret]->on_boundary )
  {
    manager.vertices_to_request().push_back( vertices[ret]->id );
    return -1;
  }
  else
    return ret;
}
//-----------------------------------------------------------------------------
bool LocalMeshCoarsening::checkMesh(std::list<DCell *>& cells_to_regenerate,
                                    std::vector<uint>& cells_to_regenerate_orient,
                                    real quality_threshold)
{
  // Check for inverted cells and new cell volumes of cells adjacent to
  // removed vertex
  std::vector<uint>::iterator o_it(cells_to_regenerate_orient.begin());
  for ( std::list<DCell *>::iterator c_it(cells_to_regenerate.begin()) ;
        c_it != cells_to_regenerate.end() ; ++c_it, ++o_it )
  {
    DCell * dc = *c_it;
    dolfin_assert( !dc->deleted );

    // check orientation of new cell
    if ( dc->orientation() != *o_it )
      return false;

    // check qm of new cell
    real qm = dc->volume() / dc->diameter();
    if ( qm <= quality_threshold )
      return false;

  }

  return true;
}
//-----------------------------------------------------------------------------
int LocalMeshCoarsening::coarsenCell(CoarseningManager& manager,
                                     DCell * cell_to_coarsen, uint attempts)
{
  // Check if cell has already been deleted
  if ( cell_to_coarsen->deleted )
    return 0;

  // Select edge for collapse
  DVertex * verts[2];
  if ( !selectEdge( cell_to_coarsen, manager, verts ) )
    return 0;

  int vert_idx = selectVertex( verts, manager, attempts );
  if ( vert_idx < 0 )
  {
    // Cannot be coarsened due to missing entities from other processes
    manager.cells_to_request().push_back( cell_to_coarsen->id );
    return -2;
  }

  DVertex * vertex_to_remove = verts[vert_idx > 0];
  DVertex * vertex_to_keep = verts[vert_idx == 0];

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
  if ( checkMesh(cells_to_regenerate, cells_to_regenerate_orientations,
                 manager.qualityThreshold()) )
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
