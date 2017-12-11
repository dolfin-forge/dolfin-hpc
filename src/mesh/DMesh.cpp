// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2013.
// Modified by Balthasar Reuter, 2013
// Modified by Aurelien Larcher, 2015
//

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/constants.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshValues.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/DMesh.h>
#include <dolfin/mesh/DVertex.h>
#include <dolfin/mesh/DCell.h>
#include <cstring>

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
DMesh::DMesh(Mesh& mesh) :
    vertices(),
    cells(),
    _is_distributed(false),
    _cell_type(NULL),
    _tdim(0),
    _gdim(0),
    _glb_max(0),
    _salt(0),
    _start_offset(0)
{
  imp(mesh);
}
//-----------------------------------------------------------------------------
DMesh::~DMesh()
{
  clear();
}
//-----------------------------------------------------------------------------
void DMesh::clear()
{
  delete _cell_type;

  // Delete allocated DCells
  for (std::list<DCell*>::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    delete *it;
  }
  cells.clear();

  // Delete allocated DVertices
  for (std::set<DVertex*>::iterator it = vertices.begin(); it != vertices.end();
      ++it)
  {
    delete *it;
  }
  vertices.clear();
}
//-----------------------------------------------------------------------------
void DMesh::init(Mesh& mesh)
{
  if (_cell_type != NULL)
  {
    error("Dynamic mesh already initialized on a mesh.");
  }

  // Cleanup before new allocation
  clear();

  _cell_type = CellType::create(mesh.type().cellType());
  _is_distributed = mesh.is_distributed();
  _tdim = mesh.topology_dimension();
  _gdim = mesh.geometry_dimension();

  // Since the mesh is linear numbered, the maximum global index assigned is
  // the number of vertices in the *global* mesh
  _glb_max = mesh.topology().global_size(0);
  dolfin_assert(_glb_max > 0);
  _salt = _cell_type->num_entities(0) * mesh.topology().global_size(_tdim);
  dolfin_assert(_salt > 0);

  // Assign a safe range for each rank for the numbering of new entities i.e
  // such that there is no overlap with existing numbered entities.
  _start_offset = 0;
  // Assume uniform refinement
  uint num_new = mesh.topology().size(1);
#ifdef HAVE_MPI
  _start_offset = 0;
#if ( MPI_VERSION > 1 )
  MPI_Exscan(&num_new, &_start_offset, 1, MPI_UNSIGNED, MPI_SUM,
             MPI::DOLFIN_COMM);
#else
  MPI_Scan(&num_new, &_start_offset, 1, MPI_UNSIGNED, MPI_SUM,
      MPI::DOLFIN_COMM);
  _start_offset -= num_new;
#endif
  _start_offset += _glb_max;
#endif
}
//-----------------------------------------------------------------------------
void DMesh::imp(Mesh& mesh)
{
  init(mesh);

  std::vector<DVertex *> vertexvec;

  // Copy vertices
  uint counter = 1;
  for (VertexIterator vi(mesh); !vi.end(); ++vi)
  {
    dolfin_assert(vi->index() == vertices.size());
    dolfin_assert(vi->index() == vertexvec.size());

    DVertex* dv = new DVertex;
    dv->p = vi->point();
    dv->id = vi->index();
    dv->glb_id = vi->global_index();
    dv->on_boundary = vi->is_shared();
    dv->shared = vi->is_shared();
    dv->ghosted = vi->is_ghost();
    if (dv->ghosted)
    {
      dv->owner = mesh.distdata()[0].get_owner(vi->index());
    }
    else if (dv->shared)
    {
      dv->shared_adj = mesh.distdata()[0].get_shared_adj(vi->index());
    }

    if (dv->on_boundary) bc_dvs[dv->glb_id] = dv;

    vertices.insert(dv);
    vertexvec.push_back(dv);
    counter++;
  }

  // Copy cells
  for (CellIterator ci(mesh); !ci.end(); ++ci)
  {
    DCell* dc = new DCell;

    std::vector<DVertex*> vs(ci->num_entities(0));
    uint i = 0;
    for (VertexIterator vi(*ci); !vi.end(); ++vi)
    {
      DVertex* dv = vertexvec[vi->index()];

      vs[i] = dv;
      i++;
    }

    add_cell(dc, vs, ci->index());
    // Define the same cell numbering
    dc->id = ci->index();

  }
}
//-----------------------------------------------------------------------------
void DMesh::exp(Mesh& mesh)
{
  eraseRemovedEntities();
  number();

  MeshEditor editor(mesh, _cell_type->cellType(), _gdim);

  editor.init_vertices(vertices.size());
  editor.init_cells(cells.size());

  // Add old vertices
  uint current_vertex = 0;
  for (std::set<DVertex*>::iterator it = vertices.begin(); it != vertices.end();
      ++it, ++current_vertex)
  {
    DVertex* dv = *it;
    dolfin_assert(!dv->deleted);

    editor.add_vertex(current_vertex, &dv->p[0]);

    if(_is_distributed)
    {
      if (dv->ghosted)
      {
        mesh.distdata()[0].set_ghost(current_vertex, dv->owner);
      }
      else if (dv->shared)
      {
        mesh.distdata()[0].set_shared(current_vertex);
      }
      mesh.distdata()[0].set_map(current_vertex, dv->glb_id);
    }
  }

  Array<uint> cell_vertices(_cell_type->num_entities(0));
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
    editor.add_cell(current_cell, &cell_vertices[0]);

    current_cell++;
  }
  editor.close();
}
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

  editor.init_vertices(vertices.size());
  editor.init_cells(cells.size());

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

    editor.add_vertex(current_vertex, &dv->p[0]);

    if(_is_distributed)
    {
      if (dv->ghosted)
      {
        mesh.distdata()[0].set_ghost(current_vertex, dv->owner);
      }
      else if (dv->shared)
      {
        mesh.distdata()[0].setall_shared_adj(current_vertex, dv->shared_adj);
      }
      mesh.distdata()[0].set_map(current_vertex, dv->glb_id);
    }
  }

  Array<uint> cell_vertices(_cell_type->num_entities(0));
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
    editor.add_cell(current_cell, &cell_vertices[0]);
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
        if (v0->glb_id > v1->glb_id) l = v0->p.dist(v1->p);
        else l = v1->p.dist(v0->p);

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
    add_vertex(mv);
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
      node.owner = MPI::rank();
      std::pair<uint, prop_edge> _prop_(dcell->nref, node);
      propagate.push_back(_prop_);
      dcell->nref++;
      bc_dvs[mv->glb_id] = mv;
      mv->on_boundary = true;
      mv->shared = true;
      mv->ghosted = false;
      mv->owner = MPI::rank();
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

  add_cell(c0, vs0, dcell->parent_id);
  add_cell(c1, vs1, dcell->parent_id);

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
void DMesh::add_vertex(DVertex* v)
{
  vertices.insert(v);
}
//-----------------------------------------------------------------------------
void DMesh::add_cell(DCell* c, std::vector<DVertex*> vs, int parent_id)
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
void DMesh::bisectMarked(MeshValues<bool, Cell> const& marked_ids)
{
  std::list<DCell*> marked_cells;
  for (std::list<DCell*>::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    DCell* c = *it;

    if (marked_ids(c->id))
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

    if (MPI::rank() == 0 && propagate.size() > 0) begin(
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

              if (MPI::rank() < it->second.owner)
              {
                mv->ghosted = false;
                mv->owner = MPI::rank();
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

    if (MPI::rank() == 0) end();

  }
}
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
  uint rank = MPI::rank();
  uint pe_size = MPI::size();
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
  uint rank = MPI::rank();
  uint pe_size = MPI::size();
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
