// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2013.
// Modified by Balthasar Reuter, 2013
//

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/constants.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/DMesh.h>
#include <dolfin/mesh/DVertex.h>
#include <dolfin/mesh/DCell.h>
#include <cstring>

#ifdef HAVE_LIBGEOM
#include <Geometry.h>
#endif

#include <algorithm>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace dolfin;
//-----------------------------------------------------------------------------
/// Helper class for getCell method
class CheckCellId
{
public:
  explicit CheckCellId(int id_) :
      id(id_)
  {
  }
  bool operator()(const DCell * const & cell) const
  {
    return (id == cell->id);
  }
private:
  int id;
};
//-----------------------------------------------------------------------------
/// Helper class for getVertex method
class CheckVertexId
{
public:
  explicit CheckVertexId(int id_) :
      id(id_)
  {
  }
  bool operator()(const DVertex * const & vertex) const
  {
    return (id == vertex->id);
  }
private:
  int id;
};
//-----------------------------------------------------------------------------
/// Helper function for eraseRemovedEntities method
static bool checkCellDeleted(DCell * cell)
{
  return cell->deleted;
}
//-----------------------------------------------------------------------------
//DMesh::DMesh() : vertices(0), cells(0)
DMesh::DMesh() :
    cells(0)
{
}
//-----------------------------------------------------------------------------
DMesh::~DMesh()
{
  if (_cell_type) delete _cell_type;

  // Delete allocated DVertices
  for (std::set<DVertex*>::iterator it = vertices.begin(); it != vertices.end();
      ++it)
    delete *it;

  // Delete allocated DCells
  for (std::list<DCell*>::iterator it = cells.begin(); it != cells.end(); ++it)
    delete *it;

}
//-----------------------------------------------------------------------------
void DMesh::imp(Mesh& mesh)
{
  _cell_type = CellType::create(mesh.type().cellType());
  //cell_type = &(mesh.type());
  _tdim = mesh.topology().dim();
  _gdim = mesh.geometry().dim();

  // Delete allocated DVertices
  for (std::set<DVertex*>::iterator it = vertices.begin(); it != vertices.end();
      ++it)
    delete *it;

  // Delete allocated DCells
  for (std::list<DCell*>::iterator it = cells.begin(); it != cells.end(); ++it)
    delete *it;

  vertices.clear();
  cells.clear();

  //MeshRenumber::renumber_vertices(mesh);

  std::vector<DVertex *> vertexvec;

  // Assume uniform refinement
  uint num_new = mesh.size(1);

  // Since the mesh is linear numbered, the maximum global index assigned is
  // the number of vertices in the mesh
  _glb_max = mesh.global_numVertices();

  Cell c(mesh, 0);
  _salt = c.numEntities(0) * mesh.global_numCells();

  // Assign a safe range for each processor
  _start_offset = 0;
#ifdef HAVE_MPI
#if ( MPI_VERSION > 1 )
  MPI_Exscan(&num_new, &_start_offset, 1, MPI_UNSIGNED, MPI_SUM,
             MPI::DOLFIN_COMM);
#else
  MPI_Scan(&num_new, &_start_offset, 1,
      MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
  _start_offset -= num_new;
#endif
  _start_offset += _glb_max;
#endif

  // Copy vertices
  uint counter = 1;
  for (VertexIterator vi(mesh); !vi.end(); ++vi)
  {
    dolfin_assert(vi->index() == vertices.size());
    dolfin_assert(vi->index() == vertexvec.size());

    DVertex* dv = new DVertex;
    dv->p = vi->point();
    dv->id = vi->index();
    dv->glb_id = mesh.distdata().get_vertex_global(vi->index());
    dv->on_boundary = MPI::numProcesses() > 1
        && mesh.distdata().is_shared(vi->index(), 0);
    dv->shared = mesh.distdata().is_shared(vi->index(), 0);
    dv->ghosted = mesh.distdata().is_ghost(vi->index(), 0);
    if (dv->ghosted) dv->owner = mesh.distdata().get_owner(*vi);
    else if (dv->shared) dv->shared_adj = mesh.distdata().get_shared_adj(*vi);

    if (dv->on_boundary) bc_dvs[dv->glb_id] = dv;

    vertices.insert(dv);
    vertexvec.push_back(dv);
    counter++;
  }

  // Copy cells
  for (CellIterator ci(mesh); !ci.end(); ++ci)
  {
    DCell* dc = new DCell;

    std::vector<DVertex*> vs(ci->numEntities(0));
    uint i = 0;
    for (VertexIterator vi(*ci); !vi.end(); ++vi)
    {
      DVertex* dv = vertexvec[vi->index()];

      vs[i] = dv;
      i++;
    }

    addCell(dc, vs, ci->index());
    // Define the same cell numbering
    dc->id = ci->index();

  }
}
//-----------------------------------------------------------------------------
#ifdef HAVE_LIBGEOM
//-----------------------------------------------------------------------------
void DMesh::imp(Mesh& mesh, MeshFunction<int>& patch_id_list,
    MeshFunction<float>& bnd_u, MeshFunction<float>& bnd_v)
{

  _cell_type = CellType::create(mesh.type().cellType());

  _tdim = mesh.topology().dim();
  _gdim = mesh.geometry().dim();

  // Delete allocated DVertices
  for(std::set<DVertex* >::iterator it = vertices.begin();
      it != vertices.end(); ++it)
  delete *it;

  // Delete allocated DCells
  for(std::list<DCell* >::iterator it = cells.begin();
      it != cells.end(); ++it)
  delete *it;

  vertices.clear();
  cells.clear();

  //MeshRenumber::renumber_vertices(mesh);

  std::vector<DVertex *> vertexvec;

  // Assume uniform refinement
  uint num_new = mesh.size(1);

  // Since the mesh is linear numbered, the maximum global index assigned is
  // the number of vertices in the mesh
  _glb_max = mesh.global_numVertices();

  Cell c(mesh, 0);
  _salt = c.numEntities(0) * mesh.distdata().global_numCells();

  // Assign a safe range for each processor
  _start_offset = 0;
#ifdef HAVE_MPI
#if ( MPI_VERSION > 1 )
  MPI_Exscan(&num_new, &_start_offset, 1,
      MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
  MPI_Scan(&num_new, &_start_offset, 1,
      MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
  _start_offset -= num_new;
#endif
  _start_offset += _glb_max;
#endif

  uint counter = 1;
  for (VertexIterator vi(mesh); !vi.end(); ++vi)
  {
    DVertex* dv = new DVertex;
    dv->p = vi->point();
    dv->glb_id = mesh.distdata().get_vertex_global(vi->index());
    dv->on_boundary = mesh.distdata().is_shared(vi->index(), 0);
    dv->shared = mesh.distdata().is_shared(vi->index(), 0);
    dv->ghosted = mesh.distdata().is_ghost(vi->index(), 0);

    dv->patch_id = patch_id_list.get(vi->index());
    dv->u = bnd_u.get(vi->index());
    dv->v = bnd_v.get(vi->index());

    if (dv->ghosted)
    dv->owner = mesh.distdata().get_owner(*vi);
    else if (dv->shared)
    dv->shared_adj = mesh.distdata().get_shared_adj(*vi);

    if(dv->on_boundary)
    bc_dvs[dv->glb_id] = dv;

    vertices.insert(dv);
    vertexvec.push_back(dv);
    counter++;
  }

  for (CellIterator ci(mesh); !ci.end(); ++ci)
  {
    DCell* dc = new DCell;

    std::vector<DVertex*> vs(ci->numEntities(0));
    uint i = 0;
    for (VertexIterator vi(*ci); !vi.end(); ++vi)
    {
      DVertex* dv = vertexvec[vi->index()];

      vs[i] = dv;
      i++;
    }

    addCell(dc, vs, ci->index());

    // Define the same cell numbering
    dc->id = ci->index();
  }
}
//-----------------------------------------------------------------------------
#endif // HAVE_LIBGEOM
//-----------------------------------------------------------------------------
void DMesh::exp(Mesh& mesh)
{
  eraseRemovedEntities();
  number();

  MeshEditor editor(mesh, _cell_type->cellType(), _gdim);

  editor.initVertices(vertices.size());
  editor.initCells(cells.size());

  // Add old vertices
  uint current_vertex = 0;
  for (std::set<DVertex*>::iterator it = vertices.begin(); it != vertices.end();
      ++it)
  {
    DVertex* dv = *it;
    dolfin_assert(!dv->deleted);

    editor.addVertex(current_vertex, dv->p);
    if (dv->ghosted)
    {
      mesh.distdata().set_ghost(current_vertex, 0);
      mesh.distdata().set_ghost_owner(current_vertex, dv->owner, 0);
    }
    else if (dv->shared) mesh.distdata().set_shared(current_vertex, 0);
    mesh.distdata().set_map(current_vertex++, dv->glb_id, 0);
  }

  Array<uint> cell_vertices(_cell_type->numEntities(0));
  uint current_cell = 0;
  for (std::list<DCell*>::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    DCell* dc = *it;
    dolfin_assert(!dc->deleted);

    for (uint j = 0; j < dc->vertices.size(); j++)
    {
      DVertex* dv = dc->vertices[j];
      cell_vertices[j] = dv->id;
    }
    editor.addCell(current_cell, cell_vertices);

    current_cell++;
  }
  editor.close();

  if (MPI::numProcesses() > 1)
  {
    mesh.distdata().set_invalid_numbering();
    mesh.distdata().set_invalid_ownership();
    mesh.renumber();
  }
}
//-----------------------------------------------------------------------------
#ifdef HAVE_LIBGEOM
//-----------------------------------------------------------------------------
void DMesh::exp(Mesh& mesh, MeshFunction<int>& patch_id_list,
    MeshFunction<float>& bnd_u, MeshFunction<float>& bnd_v)
{

  eraseRemovedEntities();
  number();

  MeshEditor editor(mesh, _cell_type->cellType(), _gdim);

  editor.initVertices(vertices.size());
  editor.initCells(cells.size());

  //initialize Meshfunctions for refined mesh
  patch_id_list.init(0);
  bnd_u.init(0);
  bnd_v.init(0);
  patch_id_list = -1;
  bnd_u = -1;
  bnd_v = -1;

  // Add old vertices
  uint current_vertex = 0;
  for(std::set<DVertex* >::iterator it = vertices.begin();
      it != vertices.end(); ++it)
  {
    DVertex* dv = *it;
    editor.addVertex(current_vertex, dv->p);
    patch_id_list.set(current_vertex, dv->patch_id);

    bnd_u.set(current_vertex, dv->u);
    bnd_v.set(current_vertex, dv->v);
    if(dv->ghosted)
    {
      mesh.distdata().set_ghost(current_vertex, 0);
      mesh.distdata().set_ghost_owner(current_vertex, dv->owner, 0);
    }
    else if(dv->shared)
    mesh.distdata().set_shared(current_vertex, 0);
    mesh.distdata().set_map(current_vertex++, dv->glb_id, 0);
  }

  Array<uint> cell_vertices(_cell_type->numEntities(0));
  uint current_cell = 0;
  for(std::list<DCell* >::iterator it = cells.begin();
      it != cells.end(); ++it)
  {
    DCell* dc = *it;

    for(uint j = 0; j < dc->vertices.size(); j++)
    {
      DVertex* dv = dc->vertices[j];
      cell_vertices[j] = dv->id;
    }
    editor.addCell(current_cell, cell_vertices);

    current_cell++;
  }
  editor.close();

  if (MPI::numProcesses() > 1)
  {
    mesh.distdata().set_invalid_numbering();
    mesh.distdata().set_invalid_ownership();
    mesh.renumber();
  }
}
//-----------------------------------------------------------------------------
#endif // HAVE_LIBGEOM
//-----------------------------------------------------------------------------
void DMesh::expKeepNumbering(Mesh& mesh, Array<int> * old2new_cells,
                             Array<int> * old2new_vertices)
{
  // Remove entities marked for deletion
  eraseRemovedEntities();

  // Renumber and create mapping
  //number(old2new_cells, old2new_vertices);

  bool delete_vertices_array(!old2new_vertices);

  //FIXME: Overallocation as the comment for glb_max was incorrect.
  if (old2new_vertices)
  dolfin_assert(old2new_vertices->size() >= vertices.size());
  else old2new_vertices = new Array<int>(_glb_max);
  *old2new_vertices = -1;

  if (old2new_cells)
  {
    dolfin_assert(old2new_cells->size() >= cells.size());
    *old2new_cells = -1;
  }

  MeshEditor editor(mesh, _cell_type->cellType(), _gdim);

  editor.initVertices(vertices.size());
  editor.initCells(cells.size());

  // Add old vertices
  uint current_vertex = 0;
  for (std::set<DVertex*>::iterator it = vertices.begin(); it != vertices.end();
      ++it, ++current_vertex)
  {
    DVertex* dv = *it;
    dolfin_assert(!dv->deleted);

#if (__sgi)
    (*old2new_vertices)[dv->id] = current_vertex;
#else
    old2new_vertices->at(dv->id) = current_vertex;
#endif

    editor.addVertex(current_vertex, dv->p);

    if (dv->ghosted)
    {
      mesh.distdata().set_ghost(current_vertex, 0);
      mesh.distdata().set_ghost_owner(current_vertex, dv->owner, 0);
    }
    else if (dv->shared)
    {
      mesh.distdata().set_shared(current_vertex, 0);
      mesh.distdata().setall_shared_adj(current_vertex, dv->shared_adj, 0);
    }

    mesh.distdata().set_map(current_vertex, dv->glb_id, 0);
  }

  Array<uint> cell_vertices(_cell_type->numEntities(0));
  uint current_cell = 0;
  for (std::list<DCell*>::iterator it = cells.begin(); it != cells.end();
      ++it, ++current_cell)
  {
    DCell* dc = *it;
    dolfin_assert(!dc->deleted);

    if (old2new_cells)
    {
#if (__sgi)
      (*old2new_cells)[dc->id] = current_cell;
#else
      old2new_cells->at(dc->id) = current_cell;
#endif
    }

    for (uint j = 0; j < dc->vertices.size(); j++)
    {
      DVertex* dv = dc->vertices[j];
#if (__sgi)
      cell_vertices[j] = (*old2new_vertices)[dv->id];
#else
      cell_vertices[j] = old2new_vertices->at(dv->id);
#endif
    }
    editor.addCell(current_cell, cell_vertices);
  }
  editor.close();

  if (delete_vertices_array) delete old2new_vertices;
}
//-----------------------------------------------------------------------------
void DMesh::number(Array<int> * old2new_cells, Array<int> * old2new_vertices)
{
  if (old2new_vertices)
  {
    dolfin_assert(old2new_vertices->size() >= vertices.size());
    *old2new_vertices = -1;
  }

  uint i = 0;
  for (std::set<DVertex*>::iterator it = vertices.begin(); it != vertices.end();
      ++it, ++i)
  {
    DVertex* dv = *it;

    if (old2new_vertices)
    {
#if (__sgi)
      (*old2new_vertices)[dv->id] = i;
#else
      old2new_vertices->at(dv->id) = i;
#endif
    }
    dv->id = i;
  }

  if (old2new_cells)
  {
    dolfin_assert(old2new_cells->size() >= cells.size());
    *old2new_cells = -1;
  }

  i = 0;
  for (std::list<DCell*>::iterator it = cells.begin(); it != cells.end();
      ++it, ++i)
  {
    DCell* dc = *it;

    if (old2new_cells)
    {
#if (__sgi)
      (*old2new_cells)[dc->id] = i;
#else
      old2new_cells->at(dc->id) = i;
#endif
    }

    dc->id = i;
  }
}
//-----------------------------------------------------------------------------
void DMesh::bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1)
{

  bool closing = false;

  // Find longest edge
  real lmax = 0.0;
  int ptmax = 0;
  uint ii = 0;
  uint jj = 0;
  for (uint i = 0; i < dcell->vertices.size(); i++)
  {
    for (uint j = 0; j < dcell->vertices.size(); j++)
    {
      if (i != j)
      {
        DVertex* v0 = dcell->vertices[i];
        DVertex* v1 = dcell->vertices[j];

        real l = 0.0;
        if (v0->glb_id > v1->glb_id) l = v0->p.distance(v1->p);
        else l = v1->p.distance(v0->p);

        if (fabs(l - lmax) < DOLFIN_EPS)
        {
          int ptsum = (v0->glb_id) + (v1->glb_id);
          if (ptsum > ptmax)
          {
            ii = i;
            jj = j;
            lmax = l;
            ptmax = (v0->glb_id + v1->glb_id);
          }
        }
        else if (l >= lmax)
        {
          ii = i;
          jj = j;
          lmax = l;
          ptmax = (v0->glb_id + v1->glb_id);
        }
      }
    }
  }

  DVertex* v0 = dcell->vertices[ii];
  DVertex* v1 = dcell->vertices[jj];
  DVertex* mv = 0;

  // Check if no hanging vertices remain, otherwise create hanging
  // vertex and continue refinement
  if ((v0 == hv0 || v0 == hv1) && (v1 == hv0 || v1 == hv1))
  {

    mv = hangv;
    closing = true;

    if (v0->on_boundary && v1->on_boundary)
    {
      mv->on_boundary = true;
      mv->shared = true;
      bc_dvs[mv->glb_id] = mv;
      dolfin_assert(v0->glb_id != v1->glb_id);
      dolfin_assert(
          ref_edge.find(edge_key(v0->glb_id, v1->glb_id)) != ref_edge.end());
    }
  }
  else
  {
    mv = new DVertex;
    addVertex(mv);
    if (v0->glb_id < v1->glb_id) mv->glb_id = (((v0->glb_id * _salt)
        + (v1->glb_id))) + _glb_max;
    else mv->glb_id = (((v1->glb_id * _salt) + (v0->glb_id))) + _glb_max;
    mv->p = (dcell->vertices[ii]->p + dcell->vertices[jj]->p) / 2.0;

    // Add hanging node on shared edges to propagation buffer
    if (v0->on_boundary && v1->on_boundary)
    {
      prop_edge node;
      node.mv = mv->glb_id;
      node.v1 = v0->glb_id;
      node.v2 = v1->glb_id;
      node.owner = MPI::processNumber();
      std::pair<uint, prop_edge> _prop_(dcell->nref, node);
      propagate.push_back(_prop_);
      dcell->nref++;
      bc_dvs[mv->glb_id] = mv;
      mv->on_boundary = true;
      mv->shared = true;
      mv->ghosted = false;
      mv->owner = MPI::processNumber();
      ref_edge[edge_key(v0->glb_id, v1->glb_id)] = mv;
    }

    closing = false;
  }

  // Create new cells
  DCell* c0 = new DCell;
  DCell* c1 = new DCell;
  c0->nref = dcell->nref;
  c1->nref = dcell->nref;
  std::vector<DVertex*> vs0(0);
  std::vector<DVertex*> vs1(0);
  for (uint i = 0; i < dcell->vertices.size(); i++)
  {
    if (i != ii) vs0.push_back(dcell->vertices[i]);

    if (i != jj) vs1.push_back(dcell->vertices[i]);
  }
  vs0.push_back(mv);
  vs1.push_back(mv);

  addCell(c0, vs0, dcell->parent_id);
  addCell(c1, vs1, dcell->parent_id);

  removeCell(dcell);

  // Continue refinement
  if (!closing)
  {
    // Bisect opposite cell of edge with hanging node
    for (;;)
    {
      DCell* copp = opposite(dcell, v0, v1);
      if (copp != 0)
      {
        bisect(copp, mv, v0, v1);
      }
      else
      {
        break;
      }
    }
  }
}
//-----------------------------------------------------------------------------
#ifdef HAVE_LIBGEOM
//-----------------------------------------------------------------------------
void DMesh::bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1,
    libgeom::Geometry& geom)
{

  bool closing = false;

  // Find longest edge
  real lmax = 0.0;
  int ptmax = 0;
  uint ii = 0;
  uint jj = 0;
  for(uint i = 0; i < dcell->vertices.size(); i++)
  {
    for(uint j = 0; j < dcell->vertices.size(); j++)
    {

      if(i != j)
      {
        DVertex* v0 = dcell->vertices[i];
        DVertex* v1 = dcell->vertices[j];

        real l = 0.0;
        if( v0->glb_id > v1->glb_id)
        l = v0->p.distance(v1->p);
        else
        l = v1->p.distance(v0->p);

        if(fabs(l - lmax) < DOLFIN_EPS)
        {
          int ptsum = (v0->glb_id) + (v1->glb_id);
          if(ptsum > ptmax)
          {
            ii = i;
            jj = j;
            lmax = l;
            ptmax = (v0->glb_id + v1->glb_id);
          }
        }
        else if(l >= lmax)
        {
          ii = i;
          jj = j;
          lmax = l;
          ptmax = (v0->glb_id + v1->glb_id);
        }
      }
    }
  }

  DVertex* v0 = dcell->vertices[ii];
  DVertex* v1 = dcell->vertices[jj];
  DVertex* mv = 0;

  // Check if no hanging vertices remain, otherwise create hanging
  // vertex and continue refinement

  if((v0 == hv0 || v0 == hv1) && (v1 == hv0 || v1 == hv1))
  {
    mv = hangv;
    closing = true;

    if( v0->on_boundary && v1->on_boundary )
    {
      mv->on_boundary = true;
      mv->shared = true;
      bc_dvs[mv->glb_id] = mv;
      dolfin_assert(v0->glb_id != v1->glb_id);
      dolfin_assert( ref_edge.find(edge_key(v0->glb_id, v1->glb_id)) != ref_edge.end());
    }
  }
  else
  {
    mv = new DVertex;
    mv->patch_id = -1;
    mv->u = -1;
    mv->v = -1;
    addVertex(mv);
    if( v0->glb_id < v1->glb_id )
    mv->glb_id = (((v0->glb_id * _salt) + (v1->glb_id))) + _glb_max;
    else
    mv->glb_id = (((v1->glb_id * _salt) + (v0->glb_id))) + _glb_max;

    //use the regular mitpoint rule when both neighbor vertices are not on the geometry
    if(dcell->vertices[ii]->patch_id < 0 || dcell->vertices[jj]->patch_id < 0)
    {
      mv->p = (dcell->vertices[ii]->p + dcell->vertices[jj]->p) / 2.0;
      mv->patch_id = -1;
      mv->u = -1;
      mv->v = -1;
    }
    else
    {
      Point midpoint = (dcell->vertices[ii]->p + dcell->vertices[jj]->p) / 2.0;
      libgeom::Point3D midpoint_lib(midpoint.x(), midpoint.y(), midpoint.z());
      libgeom::Point3D r1;
      float u1, v1;
      int pid_tmp;
      real distance;
      if(dcell->vertices[ii]->patch_id == dcell->vertices[jj]->patch_id)
      {
        pid_tmp = dcell->vertices[jj]->patch_id;
        distance = geom.find_closest_point(midpoint_lib, r1, u1, v1, pid_tmp, 100, 100, 0.00001, 0.00002, 100);
      }
      else
      {
        distance = geom.find_closest_point_all_patches(midpoint_lib, r1, u1, v1, pid_tmp, 100, 100, 8, 8, 0.00001, 0.00002, 100);
      }
      mv->p = Point( r1.x(), r1.y(), r1.z());
      mv->patch_id = pid_tmp;
      mv->u = u1;
      mv->v = v1;
    }

    // Add hanging node on shared edges to propagation buffer
    if( v0->on_boundary && v1->on_boundary)
    {
      prop_edge node;
      node.mv = mv->glb_id;
      node.v1 = v0->glb_id;
      node.v2 = v1->glb_id;
      node.owner = MPI::processNumber();
      std::pair<uint, prop_edge> _prop_( dcell->nref, node);
      propagate.push_back( _prop_ );
      dcell->nref++;
      bc_dvs[mv->glb_id] = mv;
      mv->on_boundary = true;
      mv->shared = true;
      mv->ghosted = false;
      mv->owner = MPI::processNumber();
      ref_edge[edge_key(v0->glb_id, v1->glb_id)] = mv;
    }
    //message("rank %d closing set to false", MPI::processNumber());
    closing = false;
  }

  // Create new cells
  DCell* c0 = new DCell;
  DCell* c1 = new DCell;
  c0->nref = dcell->nref;
  c1->nref = dcell->nref;
  std::vector<DVertex*> vs0(0);
  std::vector<DVertex*> vs1(0);
  for(uint i = 0; i < dcell->vertices.size(); i++)
  {
    if(i != ii)
    vs0.push_back(dcell->vertices[i]);

    if(i != jj)
    vs1.push_back(dcell->vertices[i]);
  }
  vs0.push_back(mv);
  vs1.push_back(mv);

  addCell(c0, vs0, dcell->parent_id);
  addCell(c1, vs1, dcell->parent_id);

  removeCell(dcell);

  // Continue refinement
  if(!closing)
  {
    //TODO counter stuff entfernen
    int counter = 0;
    // Bisect opposite cell of edge with hanging node
    for(;;)
    {
      counter++;
      DCell* copp = opposite(dcell, v0, v1);
      if(counter == 10000)
      {

        error("you ran into an infinite loop in the refinement. Usually this is the case if something went wrong with the initial mapping of the geometric values fromt the initial projection.");
        break;
      }
      else if(copp != 0)
      {
        bisect(copp, mv, v0, v1, geom);
      }
      else
      {
        break;
      }
    }
  }
}
//-----------------------------------------------------------------------------
#endif // HAVE_LIBGEOM
//-----------------------------------------------------------------------------
DCell* DMesh::opposite(DCell* dcell, DVertex* v1, DVertex* v2)
{
  for (std::list<DCell*>::iterator it = v1->cells.begin();
      it != v1->cells.end(); ++it)
  {
    DCell* c = *it;

    if (c != dcell && !c->deleted)
    {
      int matches = 0;
      for (uint i = 0; i < c->vertices.size(); i++)
      {
        if (c->vertices[i] == v1 || c->vertices[i] == v2)
        {
          matches++;
        }
      }

      if (matches == 2)
      {
        // Found opposite cell
        return c;
      }
    }
  }
  return 0;
}
//-----------------------------------------------------------------------------
void DMesh::addVertex(DVertex* v)
{
  vertices.insert(v);
}
//-----------------------------------------------------------------------------
void DMesh::addCell(DCell* c, std::vector<DVertex*> vs, int parent_id)
{
  for (uint i = 0; i < vs.size(); i++)
  {
    DVertex* v = vs[i];
    c->vertices.push_back(v);
    v->cells.push_back(c);
  }

  cells.push_back(c);
  c->parent_id = parent_id;
}
//-----------------------------------------------------------------------------
void DMesh::removeCell(DCell* c)
{
  c->deleted = true;
}
//-----------------------------------------------------------------------------
void DMesh::removeVertex(DVertex* v)
{
  v->deleted = true;
}
//-----------------------------------------------------------------------------
void DMesh::eraseRemovedEntities()
{
  // Remove deleted cells from global list
  for (std::list<DCell *>::iterator c_it(cells.begin()); c_it != cells.end();)
  {
    DCell * dc = *c_it;
    if (dc->deleted)
    {
      c_it = cells.erase(c_it);
      delete dc;
    }
    else ++c_it;
  }

  // Remove deleted vertices from global list
  for (std::set<DVertex *>::iterator v_it(vertices.begin());
      v_it != vertices.end(); /* blank */)
  {
    DVertex * dv = *v_it;
    if (dv->deleted)
    {
      vertices.erase(v_it++);
      delete dv;
    }
    else ++v_it;
  }
}
//-----------------------------------------------------------------------------
DVertex* DMesh::getVertex(int local_id)
{
  std::set<DVertex *>::const_iterator it = std::find_if(
      vertices.begin(), vertices.end(), CheckVertexId(local_id));

  return (it != vertices.end() ? *it : 0);
}
//-----------------------------------------------------------------------------
DCell* DMesh::getCell(int local_id)
{
  std::list<DCell *>::const_iterator it = std::find_if(cells.begin(),
                                                       cells.end(),
                                                       CheckCellId(local_id));

  return (it != cells.end() ? *it : 0);
}
//-----------------------------------------------------------------------------
void DMesh::bisectMarked(std::vector<bool> marked_ids)
{
  std::list<DCell*> marked_cells;
  for (std::list<DCell*>::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    DCell* c = *it;

    if (marked_ids[c->id])
    {
      marked_cells.push_back(c);
    }
  }

  for (std::list<DCell*>::iterator it = marked_cells.begin();
      it != marked_cells.end(); ++it)
  {
    DCell* c = *it;

    if (!c->deleted)
    {
      bisect(c, 0, 0, 0);
    }
  }

  std::vector<Propagation> propagated;
  std::list<Propagation> leftovers;

  bool empty = false;

  while (!empty)
  {

    if (MPI::processNumber() == 0 && propagate.size() > 0) begin(
        "Propagate refinement...");

    propagate_refinement(propagated, empty);

    if (empty && propagated.size() == 0) break;
    propagate.clear();

    for (std::vector<Propagation>::iterator it = propagated.begin();
        it != propagated.end(); ++it)
    {

      DVertex* mv = 0;
      dolfin_assert(it->second.v1 != it->second.v2);
      if (ref_edge.find(edge_key(it->second.v1, it->second.v2))
          != ref_edge.end())
      {
        mv = ref_edge[edge_key(it->second.v1, it->second.v2)];

        if (mv->owner > (int) it->second.owner)
        {
          mv->ghosted = true;
          mv->shared = true;
          mv->owner = it->second.owner;
        }

        continue;
      }

      DVertex* v1 = 0;
      DVertex* v2 = 0;

      if (!v1 && bc_dvs.find(it->second.v1) != bc_dvs.end())
      {
        dolfin_assert(bc_dvs.find(it->second.v1) != bc_dvs.end());
        v1 = bc_dvs[it->second.v1];
      }
      if (!v2 && bc_dvs.find(it->second.v2) != bc_dvs.end())
      {
        dolfin_assert(bc_dvs.find(it->second.v2) != bc_dvs.end());
        v2 = bc_dvs[it->second.v2];
      }

      if (!v1 || !v2)
      {
        leftovers.push_back(*it);
        continue;
      }

      for (std::list<DCell*>::iterator ic = v1->cells.begin();
          ic != v1->cells.end(); ++ic)
      {
        if (!(*ic)->deleted)
        {
          if ((*ic)->has_edge(v1, v2))
          {
            dolfin_assert((*ic)->vertices.size() > 0);
            if (mv == 0)
            {
              mv = new DVertex;
              mv->shared = true;
              mv->glb_id = it->second.mv;
              vertices.insert(mv);

              if (MPI::processNumber() < it->second.owner)
              {
                mv->ghosted = false;
                mv->owner = MPI::processNumber();
                prop_edge node;
                node.mv = mv->glb_id;
                node.v1 = it->second.v1;
                node.v2 = it->second.v2;
                node.owner = mv->owner;
                std::pair<uint, prop_edge> prop(0, node);
                propagate.push_back(prop);
              }
              else
              {
                mv->ghosted = true;
                mv->owner = it->second.owner;
              }

              mv->p = (v1->p + v2->p) / 2.0;
              mv->on_boundary = true;
              bc_dvs[mv->glb_id] = mv;
              ref_edge[edge_key(v1->glb_id, v2->glb_id)] = mv;
            }
            dolfin_assert((*ic) > 0);
            bisect((*ic), mv, v1, v2);
          }
        }
      }
    }

    propagated.clear();
    for (std::list<Propagation>::iterator it = leftovers.begin();
        it != leftovers.end(); ++it)
      propagated.push_back(*it);
    leftovers.clear();

    if (MPI::processNumber() == 0) end();

  }
}
//-----------------------------------------------------------------------------
#ifdef HAVE_LIBGEOM
//-----------------------------------------------------------------------------
void DMesh::bisectMarked(std::vector<bool> marked_ids, libgeom::Geometry& geom)
{

  std::list<DCell*> marked_cells;
  for(std::list<DCell* >::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    DCell* c = *it;

    if(marked_ids[c->id])
    {
      marked_cells.push_back(c);
    }
  }

  for(std::list<DCell* >::iterator it = marked_cells.begin();
      it != marked_cells.end(); ++it)
  {
    DCell* c = *it;

    if(!c->deleted)
    bisect(c, 0, 0, 0, geom);
  }

  std::vector<Propagation> propagated;
  std::list<Propagation> leftovers;

  bool empty = false;

  while(!empty)
  {
    if(MPI::processNumber() == 0 && propagate.size() > 0)
    begin("Propagate refinement...");

    propagate_refinement(propagated, empty);

    if( empty && propagated.size() == 0) break;
    propagate.clear();

    for(std::vector<Propagation>::iterator it = propagated.begin();
        it != propagated.end();++it)
    {
      DVertex* mv = 0;
      dolfin_assert(it->second.v1 != it->second.v2);
      if(ref_edge.find(edge_key(it->second.v1, it->second.v2)) != ref_edge.end())
      {
        mv = ref_edge[edge_key(it->second.v1, it->second.v2)];

        if( mv->owner > (int) it->second.owner)
        {
          mv->ghosted = true;
          mv->shared = true;
          mv->owner = it->second.owner;
        }
        continue;
      }

      DVertex* v1 = 0;
      DVertex* v2 = 0;

      if(!v1 && bc_dvs.find(it->second.v1) != bc_dvs.end())
      {
        dolfin_assert(bc_dvs.find(it->second.v1) != bc_dvs.end());
        v1 = bc_dvs[it->second.v1];
      }
      if(!v2 && bc_dvs.find(it->second.v2) != bc_dvs.end())
      {
        dolfin_assert(bc_dvs.find(it->second.v2) != bc_dvs.end());
        v2 = bc_dvs[it->second.v2];
      }

      if(!v1 || !v2)
      {
        leftovers.push_back(*it);
        continue;
      }

      for(std::list<DCell* >::iterator ic = v1->cells.begin();
          ic != v1->cells.end();++ic)
      {
        if(!(*ic)->deleted)
        {
          if((*ic)->has_edge(v1, v2))
          {
            dolfin_assert((*ic)->vertices.size() > 0);
            if(mv == 0)
            {
              mv = new DVertex;
              mv->shared = true;
              mv->glb_id = it->second.mv;
              vertices.insert(mv);

              if( MPI::processNumber() < it->second.owner)
              {
                mv->ghosted = false;
                mv->owner = MPI::processNumber();
                prop_edge node;
                node.mv = mv->glb_id;
                node.v1 = it->second.v1;
                node.v2 = it->second.v2;
                node.owner = mv->owner;
                std::pair<uint, prop_edge> prop(0, node);
                propagate.push_back(prop);
              }
              else
              {
                mv->ghosted = true;
                mv->owner = it->second.owner;
              }

              //use the regular midpoint rule when both neighbor vertices are not on the geometry
              if(v1->patch_id < 0 || v2->patch_id < 0)
              {

                mv->p = (v1->p + v2->p) / 2.0;
                mv->patch_id = -1;
                mv->u = -1;
                mv->v = -1;
              }
              else
              {
                Point midpoint = (v1->p + v2->p) / 2.0;
                libgeom::Point3D midpoint_lib(midpoint.x(), midpoint.y(), midpoint.z());
                libgeom::Point3D r1;
                float um, vm;
                int pid_tmp;
                real distance;

                if(v1->patch_id == v2->patch_id)
                {
                  pid_tmp = v1->patch_id;
                  distance = geom.find_closest_point(midpoint_lib, r1, um, vm, pid_tmp, 100, 100, 0.00001, 0.00002, 100);
                }
                else
                {
                  distance = geom.find_closest_point_all_patches(midpoint_lib, r1, um, vm, pid_tmp, 100, 100, 8, 8, 0.00001, 0.00002, 100);
                }
                mv->p = Point( r1.x(), r1.y(), r1.z());
                mv->patch_id = pid_tmp;
                mv->u = um;
                mv->v = vm;

              }
              mv->on_boundary = true;
              bc_dvs[mv->glb_id] = mv;
              ref_edge[edge_key(v1->glb_id, v2->glb_id)] = mv;
            }
            dolfin_assert((*ic) > 0);
            bisect((*ic), mv, v1, v2, geom);
          }
        }
      }
    }

    propagated.clear();
    for(std::list<Propagation>::iterator it = leftovers.begin(); it != leftovers.end(); ++it)
    propagated.push_back(*it);
    leftovers.clear();

    if(MPI::processNumber() == 0)
    end();

  }
}
//-----------------------------------------------------------------------------
#endif // HAVE_LIBGEOM
//-----------------------------------------------------------------------------
#ifdef HAVE_MPI
//-----------------------------------------------------------------------------
void DMesh::propagate_naive(std::vector<Propagation>& propagated, bool& empty)
{
  // Allocate receive buffer
  int num_prop = propagate.size() * 5;
  int max_prop, recv_count;
  MPI_Allreduce(&num_prop, &max_prop, 1, MPI_INTEGER, MPI_MAX,
                MPI::DOLFIN_COMM);

  int *recv_buff = new int[max_prop];
  int *send_buff = new int[num_prop];
  int *sp = &send_buff[0];

  for (std::vector<Propagation>::iterator it = propagate.begin();
      it != propagate.end(); ++it)
  {
    *(sp++) = it->first;
    *(sp++) = it->second.mv;
    *(sp++) = it->second.v1;
    *(sp++) = it->second.v2;
    *(sp++) = it->second.owner;
  }

  MPI_Status status;
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint dest, src;

  empty = true;
  for (uint j = 1; j < pe_size; j++)
  {
    src = (rank - j + pe_size) % pe_size;
    dest = (rank + j) % pe_size;

    MPI_Sendrecv(&send_buff[0], num_prop, MPI_INTEGER, dest, 1, recv_buff,
                 max_prop, MPI_INTEGER, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_INTEGER, &recv_count);

    if (recv_count > 0) empty = false;

    dolfin_assert(recv_count % 5 == 0);
    for (int k = 0; k < recv_count; k += 5)
    {

      prop_edge node;
      node.mv = recv_buff[k + 1];
      node.v1 = recv_buff[k + 2];
      node.v2 = recv_buff[k + 3];
      node.owner = recv_buff[k + 4];

      Propagation prop(recv_buff[k], node);
      propagated.push_back(prop);
    }

  }

  less_pair comp;
  std::sort(propagated.begin(), propagated.end(), comp);

  short prop, gprop;
  prop = (empty == false);
  MPI_Allreduce(&prop, &gprop, 1, MPI_SHORT, MPI_SUM, MPI::DOLFIN_COMM);
  empty = (gprop == 0);

  delete[] send_buff;
  delete[] recv_buff;
}
//-----------------------------------------------------------------------------
void DMesh::propagate_hypercube(std::vector<Propagation>& propagated,
                                bool& empty)
{

  // Allocate receive buffer
  int num_prop = propagate.size() * 5;
  int total_prop, recv_count;
  MPI_Allreduce(&num_prop, &total_prop, 1, MPI_INTEGER, MPI_SUM,
                MPI::DOLFIN_COMM);

  int *recv_buff = new int[total_prop];
  int *state = new int[total_prop];
  int *sp = &state[0];
  uint state_size = 0;

  for (std::vector<Propagation>::iterator it = propagate.begin();
      it != propagate.end(); ++it)
  {
    *(sp++) = it->first;
    *(sp++) = it->second.mv;
    *(sp++) = it->second.v1;
    *(sp++) = it->second.v2;
    *(sp++) = it->second.owner;
    state_size += 5;
  }

  MPI_Status status;
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint dest;
  uint D = 1;
#if  (__sgi || __FreeBSD__)
  uint _log2, x;
  x = pe_size;
  _log2 = 0;
  while(x > 1)
  { _log2++; x>>=1;}
  for(uint j = 0; j < _log2; j++)
#else
  for (uint j = 0; j < log2(pe_size); j++)
#endif
  {
    dest = rank ^ (D << j);

    MPI_Sendrecv(state, state_size, MPI_INTEGER, dest, 1, recv_buff, total_prop,
                 MPI_INTEGER, dest, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_INTEGER, &recv_count);

    dolfin_assert(recv_count % 5 == 0);
    for (int k = 0; k < recv_count; k += 5)
    {

      prop_edge node;
      node.mv = recv_buff[k + 1];
      node.v1 = recv_buff[k + 2];
      node.v2 = recv_buff[k + 3];
      node.owner = recv_buff[k + 4];

      Propagation prop(recv_buff[k], node);
      propagated.push_back(prop);
    }
    memcpy(sp, recv_buff, recv_count * sizeof(int));
    sp += recv_count;
    state_size += recv_count;

  }

  less_pair comp;
  std::sort(propagated.begin(), propagated.end(), comp);
  empty = (state_size == 0);

  delete[] recv_buff;
  delete[] state;
}
//-----------------------------------------------------------------------------
#else
//-----------------------------------------------------------------------------
void DMesh::propagate_naive(std::vector<Propagation>& propagated, bool& empty)
{
  error("Rivara needs MPI");
}
//-----------------------------------------------------------------------------
void DMesh::propagate_hypercube(std::vector<Propagation>& propagated, bool& empty)
{
  error("Rivara needs MPI");
}
//-----------------------------------------------------------------------------
#endif
