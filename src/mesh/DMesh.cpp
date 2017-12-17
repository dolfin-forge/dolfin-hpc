// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2013.
// Modified by Balthasar Reuter, 2013
// Modified by Aurelien Larcher, 2015
//

#include <dolfin/mesh/DMesh.h>

#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/MeshValues.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/DVertex.h>
#include <dolfin/mesh/DCell.h>

#include <algorithm>
#include <cstring>

using namespace dolfin;
//-----------------------------------------------------------------------------
/// Helper class
struct CheckId
{
  explicit CheckId(int id_) : id(id_) {}
  bool operator()(const DVertex * const & v) const { return (id == v->id); }
  bool operator()(const DCell * const & c) const   { return (id == c->id); }

private:

  int const id;
};
//-----------------------------------------------------------------------------
DMesh::DMesh(Mesh& mesh) :
    vertices(),
    cells(),
    mesh_(mesh),
    ctype_(mesh.type().clone()),
    space_(mesh.space().clone()),
    glb_max_(mesh.global_size(0)),
    salt_(ctype_->num_entities(0) * mesh.global_size(mesh.type().dim())),
    start_offset_(0)
{
  dolfin_assert(glb_max_ > 0);
  dolfin_assert(salt_ > 0);

  // Assign a safe range for each rank for the numbering of new entities i.e
  // such that there is no overlap with existing numbered entities.
#ifdef HAVE_MPI
  MPI::offset(mesh.size(1), start_offset_, MPI::DOLFIN_COMM);
#endif
  start_offset_ += glb_max_;

  DVertex ** vertices = (mesh.size(0) ? new DVertex *[mesh.size(0)] : NULL);

  // Copy vertices
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    DVertex* dv = new DVertex(*v);

    if (dv->on_boundary) bc_dvs[dv->glb_id] = dv;

    add_vertex(dv);

    vertices[v->index()] = dv;
  }

  // Copy cells
  std::vector<DVertex*> vs(mesh.type().num_entities(0));
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    DCell* dc = new DCell(*c);
    uint i = 0;
    for (VertexIterator vi(*c); !vi.end(); ++vi, ++i)
    {
      vs[i] = vertices[vi->index()];
    }
    add_cell(dc, vs, c->index());
  }

  delete[] vertices;
}
//-----------------------------------------------------------------------------
DMesh::~DMesh()
{
  delete ctype_;
  delete space_;

  // Delete allocated DCells
  for (CellList::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    delete *it;
  }
  // Delete allocated DVertices
  for (VertexSet::iterator it = vertices.begin(); it != vertices.end(); ++it)
  {
    delete *it;
  }
}
//-----------------------------------------------------------------------------
void DMesh::exp(Mesh& mesh)
{
  eraseRemovedEntities();
  number();

  MeshEditor editor(mesh, *ctype_, *space_);

  editor.init_vertices(vertices.size());
  editor.init_cells(cells.size());

  // Add old vertices
  DistributedData * const dist = (mesh.is_distributed() ? &mesh.distdata()[0]
                                                        : NULL);
  uint current_vertex = 0;
  for (VertexSet::iterator it = vertices.begin(); it != vertices.end(); ++it,
       ++current_vertex)
  {
    DVertex * const dv(*it);
    dolfin_assert(!dv->deleted);

    editor.add_vertex(current_vertex, &dv->p[0]);

    if (dist)
    {
      if (dv->shared)
      {
        dist->setall_shared_adj(current_vertex, dv->shared_adj);
        if (dv->ghosted)
        {
          dist->set_ghost(current_vertex, dv->owner);
        }
      }
      dist->set_map(current_vertex, dv->glb_id);
    }
  }

  Array<uint> cell_vertices(ctype_->num_entities(0));
  uint current_cell = 0;
  for (CellList::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    DCell * const dc(*it);
    dolfin_assert(!dc->deleted);

    for (uint j = 0; j < dc->vertices.size(); j++)
    {
      DVertex* dv = dc->vertices[j];
      cell_vertices[j] = dc->vertices[j]->id;
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
  {
    dolfin_assert(old2new_vertices->size() >= vertices.size());
  }
  else
  {
    warning("Fix the allocation: the size is the global number of vertices!");
    old2new_vertices = new Array<int>(glb_max_);
  }
  *old2new_vertices = -1;

  if (old2new_cells)
  {
    dolfin_assert(old2new_cells->size() >= cells.size());
    *old2new_cells = -1;
  }

  MeshEditor editor(mesh, *ctype_, *space_);

  editor.init_vertices(vertices.size());
  editor.init_cells(cells.size());

  // Add old vertices
  DistributedData * const dist = (mesh.is_distributed() ? &mesh.distdata()[0]
                                                        : NULL);
  uint current_vertex = 0;
  for (VertexSet::iterator it = vertices.begin(); it != vertices.end(); ++it,
       ++current_vertex)
  {
    DVertex* dv = *it;
    dolfin_assert(!dv->deleted);

    (*old2new_vertices)[dv->id] = current_vertex;

    editor.add_vertex(current_vertex, &dv->p[0]);

    if (dist)
    {
      if (dv->shared)
      {
        dist->setall_shared_adj(current_vertex, dv->shared_adj);
        if (dv->ghosted)
        {
          dist->set_ghost(current_vertex, dv->owner);
        }
      }
      dist->set_map(current_vertex, dv->glb_id);
    }
  }

  Array<uint> cell_vertices(ctype_->num_entities(0));
  uint current_cell = 0;
  for (CellList::iterator it = cells.begin(); it != cells.end();
      ++it, ++current_cell)
  {
    DCell* dc = *it;
    dolfin_assert(!dc->deleted);

    if (old2new_cells)
    {
      (*old2new_cells)[dc->id] = current_cell;
    }

    for (uint j = 0; j < dc->vertices.size(); j++)
    {
      DVertex* dv = dc->vertices[j];
      cell_vertices[j] = (*old2new_vertices)[dv->id];
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
  for (VertexSet::iterator it = vertices.begin(); it != vertices.end(); ++it,
       ++i)
  {
    DVertex* dv = *it;

    if (old2new_vertices)
    {
      (*old2new_vertices)[dv->id] = i;
    }
    dv->id = i;
  }

  if (old2new_cells)
  {
    dolfin_assert(old2new_cells->size() >= cells.size());
    *old2new_cells = -1;
  }

  i = 0;
  for (CellList::iterator it = cells.begin(); it != cells.end();
       ++it, ++i)
  {
    DCell* dc = *it;

    if (old2new_cells)
    {
      (*old2new_cells)[dc->id] = i;
    }

    dc->id = i;
  }
}
//-----------------------------------------------------------------------------
void DMesh::bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1)
{
  uint const pe_rank = PE::rank();
  bool closing = false;

  // Find longest edge
  real lmax = 0.0;
  int ptmax = 0;
  uint ii = 0;
  uint jj = 0;
  for (uint i = 0; i < dcell->vertices.size(); ++i)
  {
    DVertex* const v0 = dcell->vertices[i];
    for (uint j = i + 1; j < dcell->vertices.size(); +j)
    {
      DVertex* const v1 = dcell->vertices[j];

      real const l = v0->p.dist(v1->p);
      int const ptsum = (v0->glb_id) + (v1->glb_id);
      if ((l >= lmax + DOLFIN_EPS) ||
          ((ptsum > ptmax) && fabs(l - lmax) < DOLFIN_EPS))
      {
        ii = i;
        jj = j;
        lmax = l;
        ptmax = ptsum;
      }
    }
  }

  DVertex* const v0 = dcell->vertices[ii];
  DVertex* const v1 = dcell->vertices[jj];
  DVertex* mv = NULL;

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
      //Fix shared_adj
      mv->shared_adj = ;
      bc_dvs[mv->glb_id] = mv;
      dolfin_assert(v0->glb_id != v1->glb_id);
      dolfin_assert(ref_edge.find(edge_key(v0->glb_id, v1->glb_id)) != ref_edge.end());
    }
  }
  else
  {
    mv = new DVertex;
    add_vertex(mv);
    if (v0->glb_id < v1->glb_id)
    {
      mv->glb_id = (((v0->glb_id * salt_) + (v1->glb_id))) + glb_max_;
    }
    else
    {
      mv->glb_id = (((v1->glb_id * salt_) + (v0->glb_id))) + glb_max_;
    }
    mv->p = (dcell->vertices[ii]->p + dcell->vertices[jj]->p) / 2.0;

    // Add hanging node on shared edges to propagation buffer
    if (v0->on_boundary && v1->on_boundary)
    {
      prop_edge node;
      node.mv = mv->glb_id;
      node.v1 = v0->glb_id;
      node.v2 = v1->glb_id;
      node.owner = pe_rank;
      std::pair<uint, prop_edge> _prop_(dcell->nref, node);
      propagate.push_back(_prop_);
      dcell->nref++;
      bc_dvs[mv->glb_id] = mv;
      mv->on_boundary = true;
      mv->shared = true;
      //Fix shared_adj
      mv->shared_adj = ;
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
  for (CellList::iterator it = v1->cells.begin();
      it != v1->cells.end(); ++it)
  {
    DCell* const c = *it;

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
  for (std::list<DCell *>::iterator it(cells.begin()); it != cells.end();)
  {
    if ((*it)->deleted) { delete (*it); it = cells.erase(it); } else ++it;
  }

  // Remove deleted vertices from global list
  for (std::set<DVertex *>::iterator it(vertices.begin());
       it != vertices.end(); ++it)
  {
    if ((*it)->deleted) { delete (*it); vertices.erase(it); }
  }
}
//-----------------------------------------------------------------------------
DVertex* DMesh::getVertex(int local_id)
{
  std::set<DVertex *>::const_iterator it =
      std::find_if(vertices.begin(), vertices.end(), CheckId(local_id));

  return (it != vertices.end() ? *it : NULL);
}
//-----------------------------------------------------------------------------
DCell* DMesh::getCell(int local_id)
{
  std::list<DCell *>::const_iterator it =
      std::find_if(cells.begin(), cells.end(), CheckId(local_id));

  return (it != cells.end() ? *it : NULL);
}
//-----------------------------------------------------------------------------
void DMesh::bisectMarked(MeshValues<bool, Cell> const& marked_ids)
{
  uint const pe_rank = PE::rank();

  for (CellList::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    if (marked_ids((*it)->id) && !(*it)->deleted)
    {
      bisect((*it), NULL, NULL, NULL);
    }
  }

  std::vector<Propagation> propagated;
  std::list<Propagation> leftovers;

  bool empty = false;

  while (!empty)
  {

    if (pe_rank == 0 && propagate.size() > 0)
    {
      begin("Propagate refinement...");
    }

    propagate_refinement(propagated, empty);

    if (empty && propagated.size() == 0) break;
    propagate.clear();

    for (std::vector<Propagation>::iterator it = propagated.begin();
        it != propagated.end(); ++it)
    {

      DVertex* mv = NULL;
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

      DVertex* v1 = NULL;
      DVertex* v2 = NULL;

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

      for (CellList::iterator ic = v1->cells.begin(); ic != v1->cells.end();
           ++ic)
      {
        if (!(*ic)->deleted && (*ic)->has_edge(v1, v2))
        {
          dolfin_assert((*ic)->vertices.size() > 0);
          if (mv == NULL)
          {
            mv = new DVertex;
            mv->shared = true;
            mv->glb_id = it->second.mv;
            vertices.insert(mv);

            if (pe_rank < it->second.owner)
            {
              mv->ghosted = false;
              mv->owner = pe_rank;
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
            dolfin_assert(v1->glb_id != v2->glb_id);
            ref_edge[edge_key(v1->glb_id, v2->glb_id)] = mv;
          }
          dolfin_assert((*ic) > 0);
          bisect((*ic), mv, v1, v2);
        }
      }
    }

    propagated.assign(leftovers.begin(), leftovers.end());
    leftovers.clear();

    if (pe_rank == 0) end();

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
    MPI_INTEGER,
                 dest, 1, MPI::DOLFIN_COMM, &status);
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
