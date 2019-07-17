// Copyright (C) 2008 Johan Jansson
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2009-2013, 2018.
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

#include <algorithm>
#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
struct DCell
{

  DCell() :
      id(0),
      parent_id(0),
      vertices(0),
      deleted(false),
      nref(0)
  {
  }

  DCell(Cell const& c) :
      id(c.index()),
      parent_id(0),
      vertices(0),
      deleted(false),
      nref(0)
  {
  }

  bool has_edge(DVertex *v1, DVertex *v2)
  {
    uint found = 0;
    for (std::vector<DVertex*>::iterator it = vertices.begin();
         it != vertices.end(); ++it )
      if ( (*it == v1 || *it == v2) && (++found == 2)) return true;
    return false;
  }

  /// Local index of cell
  int id;

  /// Index of parent cell
  int parent_id;

  /// List of vertices spaning the cell
  std::vector<DVertex *> vertices;

  /// Marker for deletion
  bool deleted;

  /// reference number to be used for identification in bisect
  int nref;
};
//-----------------------------------------------------------------------------
struct DVertex
{
  static uint const UNDEF = DOLFIN_UINT_MAX;

  DVertex() :
      id(UNDEF),
      glb_id(DOLFIN_LONG_MAX),
      cells(),
      p(),
      deleted(false),
      is_shared(false),
      owner(UNDEF)
  {
  }

  DVertex(Vertex const& v) :
      id(v.index()),
      glb_id(static_cast<long>(v.global_index())),
      cells(0),
      p(v.point()),
      deleted(false),
      is_shared(v.is_shared()),
      owner(v.owner())
  {
  }

  /// Local index of vertex
  uint id;

  /// Global index of vertex
  long glb_id;

  /// List of cells containing the vertex
  std::list<DCell *> cells;

  /// Vertex coordinates as Point object
  Point p;

  /// Marker for deletion
  bool deleted;

  /// Marker for shared vertices
  bool is_shared;

  /// Rank of owning process
  uint owner;

};
//-----------------------------------------------------------------------------
void sanitize_check(Mesh& mesh)
{
  // Check distributed data consistency
  if (mesh.is_distributed())
  {
    mesh.distdata()[0].check_ghost();
    mesh.distdata()[0].check_shared();
  }
}
//-----------------------------------------------------------------------------
DMesh::DMesh(Mesh& mesh) :
    vertices(),
    cells(),
    mesh_(mesh),
    ctype_(mesh.type().clone()),
    space_(mesh.space().clone()),
    shared_edges_(mesh.is_distributed() ? new SharedEdges() : NULL),
    glb_max_(static_cast<long>(mesh.global_size(0))),
    salt_((static_cast<long>(ctype_->num_entities(0)) * 
	   static_cast<long>(mesh.global_size(mesh.type().dim())))),
    cdeleted_(0),
    vdeleted_(0)
{
  dolfin_assert(glb_max_ > 0);
  dolfin_assert(salt_ > 0);
#if DEBUG
  message("DMesh: import sanitize check");
  sanitize_check(mesh);
#endif

  DVertex ** vertices = (mesh.size(0) ? new DVertex *[mesh.size(0)] : NULL);

  // Copy vertices
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    DVertex* dv = new DVertex(*v);

    if (dv->is_shared) bc_dvs[dv->glb_id] = dv;

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

  // Make sure edges are created
  //  mesh.init(1);
  mesh.init();

  // Cache shared edges
  if (shared_edges_)
  {
    DistributedData const& vdist = mesh.distdata()[0];
#if DEBUG
    {
      message("DMesh: shared edge creation sanitize check");
      mesh.distdata()[1].check_ghost();
      mesh.distdata()[1].check_shared();
    }
#endif
    for (Edge::shared it(mesh); it.valid(); ++it)
    {
      dolfin_assert(!it.adj().empty());
      Edge e(mesh, it.index());
      uint const * v = e.entities(0);
      EdgeKey<long> key(static_cast<long>(vdist.get_global(v[0])), 
			static_cast<long>(vdist.get_global(v[1])));
      dolfin_assert(vdist.is_shared(v[0]));
      dolfin_assert(vdist.is_shared(v[1]));
      shared_edges_->insert(SharedEdgeItem(key, e.index()));
    }
  }

  delete[] vertices;
}
//-----------------------------------------------------------------------------
DMesh::~DMesh()
{
  delete shared_edges_;
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
  _map<long, uint> new_global;
  {
    // Remove deleted cells from global list
    if (cdeleted_)
    {
      for (CellList::iterator it(cells.begin()); it != cells.end();)
      {
        if ((*it)->deleted) { delete (*it); it = cells.erase(it); } else ++it;
      }
    }
    // Renumber
    uint i = 0;
    for (CellList::iterator it = cells.begin(); it != cells.end();
         ++it, ++i)
    {
      (*it)->id = i;
    }

    // Renumber glb indicies    
    renumber_glb(new_global);
  }
  {
    // Remove deleted vertices from global list
    if (vdeleted_)
    {
      for (std::set<DVertex *>::iterator it(vertices.begin());
           it != vertices.end(); ++it)
      {
        if ((*it)->deleted) { delete (*it); vertices.erase(it); }
      }
    }
    // Renumber
    uint i = 0;
    for (VertexSet::iterator it = vertices.begin(); it != vertices.end(); ++it,
         ++i)
    {
      (*it)->id = i;
    }
  }

  MeshEditor editor(mesh, *ctype_, *space_);

  editor.init_vertices(vertices.size());
  editor.init_cells(cells.size());

  // Add old vertices
  DistributedData * const dist = (shared_edges_ ? &mesh.distdata()[0] : NULL);
  uint const pe_rank = mesh.topology().comm_rank();
  uint current_vertex = 0;
  for (VertexSet::iterator it = vertices.begin(); it != vertices.end(); ++it,
       ++current_vertex)
  {
    DVertex * const dv(*it);
    dolfin_assert(!dv->deleted);

    editor.add_vertex(current_vertex, &dv->p[0]);

    if (dist)
    {
      dist->set_map(current_vertex, new_global[dv->glb_id]);
      if (dv->is_shared)
      {
        dolfin_assert(dv->owner != DVertex::UNDEF);
	
	if (dv->owner != pe_rank)
	{
	  dist->set_ghost(current_vertex, dv->owner);	  
	}

      }      
    }
  }
  dist->remap_shared_adj();
  dist->finalize();

  Array<uint> cell_vertices(ctype_->num_entities(0));
  uint current_cell = 0;
  for (CellList::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    DCell * const dc(*it);
    dolfin_assert(!dc->deleted);

    for (uint j = 0; j < dc->vertices.size(); j++)
    {
      cell_vertices[j] = dc->vertices[j]->id;
    }
    editor.add_cell(current_cell, &cell_vertices[0]);

    current_cell++;
  }
  editor.close();

#if DEBUG
  message("DMesh: export sanitize check");
  sanitize_check(mesh);
#endif

}
//-----------------------------------------------------------------------------
void DMesh::bisect(DCell* dcell, DVertex* hangv, DVertex* hv0, DVertex* hv1)
{
  bool closing = false;

  // Find longest edge
  real lmax = 0.0;
  long ptmax = 0;
  DVertex * v0 = NULL;
  DVertex * v1 = NULL;
  DVertex ** const vb = &dcell->vertices[0];
  DVertex ** const ve = vb + dcell->vertices.size();
  for (DVertex ** vi = vb; vi != ve; ++vi)
  {
    for (DVertex ** vj = vi + 1; vj != ve; ++vj)
    {
      real l = 0.0;
      if ( (*vi)->glb_id > (*vj)->glb_id)
	l = (*vi)->p.dist((*vj)->p);
      else
	l = (*vj)->p.dist((*vi)->p);

      long ptsum = ((*vi)->glb_id + (*vj)->glb_id);
      if (fabs(l - lmax) < DOLFIN_EPS)
      {
	if (ptsum > ptmax) 
	{
	  v0 = *vi;
	  v1 = *vj;
	  lmax = l;
	  ptmax = ptsum;
	}
      }
      else if (l >= lmax) {
	v0 = *vi;
	v1 = *vj;
	lmax = l;
	ptmax = ptsum;
      }
    }
  }
  dolfin_assert(v0 != v1);

  DVertex* mv = NULL;

  // Check if no hanging vertices remain, otherwise create hanging
  // vertex and continue refinement
  if ((v0 == hv0 || v0 == hv1) && (v1 == hv0 || v1 == hv1))
  {

    mv = hangv;
    closing = true;

    // Having both vertices on the boundary is not a sufficient condition to be
    // shared so we need to lookup in the set of shared edges.
    if (v0->is_shared && v1->is_shared)
    {
      dolfin_assert(v0->glb_id != v1->glb_id);
      EdgeKey<long> key(v0->glb_id, v1->glb_id);
      /// @todo re-enable this check
      //SharedEdges::const_iterator it = shared_edges_->find(key);
      //      if (it != shared_edges_->end())
      //      {
	bc_dvs[mv->glb_id] = mv;
	mv->is_shared = true;
	dolfin_assert(mv->owner != DVertex::UNDEF);
	dolfin_assert(ref_edge.find(key) != ref_edge.end());
	//      }
    }
  }
  else
  {
    mv = new DVertex();
    add_vertex(mv);
    if (v0->glb_id < v1->glb_id)
    {
      dolfin_assert((v0->glb_id * salt_)
                      < (DOLFIN_LONG_MAX - glb_max_ - v1->glb_id));
      mv->glb_id = (((v0->glb_id * salt_) + (v1->glb_id))) + glb_max_;
    }
    else
    {
      dolfin_assert((v1->glb_id * salt_)
                      < (DOLFIN_LONG_MAX - glb_max_ - v0->glb_id));
      mv->glb_id = (((v1->glb_id * salt_) + (v0->glb_id))) + glb_max_;
    }
    dolfin_assert(mv->glb_id > glb_max_);
    mv->p = (v0->p + v1->p) / 2.0;

    // Add hanging node on shared edges to propagation buffer
    // Unfortunately this is a necessary condition but not sufficient.
    if (v0->is_shared && v1->is_shared)
    {
      EdgeKey<long> key(v0->glb_id, v1->glb_id);
      /// @todo re-enable this check and create new shared edges
      //      SharedEdges::const_iterator it = shared_edges_->find(key);
      //if (it != shared_edges_->end())
      //{
	prop_edge node;
	node.mv = mv->glb_id;
	node.v1 = v0->glb_id;
	node.v2 = v1->glb_id;
	node.owner = mesh_.topology().comm_rank();
	std::pair<uint, prop_edge> _prop_(dcell->nref, node);
	propagate.push_back(_prop_);
	dcell->nref++;
	bc_dvs[mv->glb_id] = mv;
	mv->is_shared = true;
	mv->owner = node.owner;
	ref_edge[key] = mv;
	//      }
    }
    closing = false;
  }

  // Create new cells
  DCell* c0 = new DCell();
  DCell* c1 = new DCell();
  c0->nref = dcell->nref;
  c1->nref = dcell->nref;
  std::vector<DVertex*> vs0;
  std::vector<DVertex*> vs1;
  for (DVertex ** vi = vb; vi != ve; ++vi)
  {
    if (*vi != v0) vs0.push_back(*vi);
    if (*vi != v1) vs1.push_back(*vi);
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
      if (copp != NULL)
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
  for (std::list<DCell *>::iterator c = v1->cells.begin(); c != v1->cells.end();
       ++c)
  {
    if ((*c) == dcell || (*c)->deleted) continue;

    for (std::vector<DVertex *>::const_iterator vi = (*c)->vertices.begin();
         vi != (*c)->vertices.end(); ++vi)
    {
      if (*vi == v2) return *c;
    }
  }
  return NULL;
}
//-----------------------------------------------------------------------------
void DMesh::add_vertex(DVertex* v)
{
  vertices.insert(v);
}
//-----------------------------------------------------------------------------
void DMesh::add_cell(DCell* c, std::vector<DVertex*> vs, int parent_id)
{
  c->vertices.assign(vs.begin(), vs.end());
  for (std::vector<DVertex*>::iterator it = vs.begin(); it != vs.end(); ++it)
  {
    (*it)->cells.push_back(c);
  }
  c->id = cells.size();
  cells.push_back(c);
  c->parent_id = parent_id;
}
//-----------------------------------------------------------------------------
void DMesh::removeCell(DCell* c)
{
  c->deleted = true;
  ++cdeleted_;
}
//-----------------------------------------------------------------------------
void DMesh::removeVertex(DVertex* v)
{
  v->deleted = true;
  ++vdeleted_;
}
//-----------------------------------------------------------------------------
void DMesh::bisectMarked(MeshValues<bool, Cell> const& marked_ids)
{
  uint const pe_rank = mesh_.topology().comm_rank();

  uint const numcells = cells.size();
  for (CellList::iterator it = cells.begin(); it != cells.end(); ++it)
  {
    if ((*it)->id >= numcells) break;
    if (!(*it)->deleted && (marked_ids((*it)->id)))
    {
      bisect((*it), NULL, NULL, NULL);
    }
  }

  Array<Propagation> propagated;
  std::list<Propagation> leftovers;

  bool empty = !marked_ids.mesh().is_distributed();

  begin("Rivara refinement");
  while (!empty)
  {

    if (pe_rank == 0 && propagate.size() > 0)
    {
      message("Propagate refinement...");
    }

    propagate_refinement(marked_ids.mesh(), propagated, empty);

    if (empty && propagated.size() == 0) break;
    propagate.clear();

    for (Array<Propagation>::iterator it = propagated.begin();
        it != propagated.end(); ++it)
    {

      DVertex* mv = NULL;
      dolfin_assert(it->second.v1 != it->second.v2);
      EdgeKey<long> key(it->second.v1, it->second.v2);
      RefinedEdges::iterator re = ref_edge.find(key);
      if (re != ref_edge.end())
      {
        mv = re->second;
	dolfin_assert(mv->is_shared);

	if (mv->owner > static_cast<uint>(it->second.owner))
	{
	  mv->owner = it->second.owner;
	  mv->is_shared = true;
	  dolfin_assert(mv->owner != DVertex::UNDEF);
	}
	
        continue;
      }

      BoundaryVertices::iterator v1it = bc_dvs.find(it->second.v1);
      if (v1it == bc_dvs.end()) { leftovers.push_back(*it); continue; }
      BoundaryVertices::iterator v2it = bc_dvs.find(it->second.v2);
      if (v2it == bc_dvs.end()) { leftovers.push_back(*it); continue; }

      DVertex* const v1 = v1it->second;
      DVertex* const v2 = v2it->second;
      dolfin_assert(v1->is_shared && v2->is_shared);

      for (std::list<DCell *>::iterator ic = v1->cells.begin(); ic != v1->cells.end();
           ++ic)
      {
        if (!(*ic)->deleted && (*ic)->has_edge(v1, v2))
        {
          dolfin_assert((*ic)->vertices.size() > 0);
          if (mv == NULL)
          {
            mv = new DVertex();
            mv->glb_id = it->second.mv;
	    mv->is_shared = true;
            vertices.insert(mv);

            if (pe_rank < it->second.owner)
            {
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
              mv->owner = it->second.owner;
            }

            mv->p = (v1->p + v2->p) / 2.0;
            bc_dvs[mv->glb_id] = mv;
            dolfin_assert(v1->glb_id != v2->glb_id);
            ref_edge[EdgeKey<long>(v1->glb_id, v2->glb_id)] = mv;
          }
          dolfin_assert((*ic) > 0);
          bisect((*ic), mv, v1, v2);
        }
      }
    }

    propagated.assign(leftovers.begin(), leftovers.end());
    leftovers.clear();



  }
  if (pe_rank == 0) end();
}

//-----------------------------------------------------------------------------
void DMesh::propagate_refinement(Mesh& mesh, Array<Propagation>& propagated, bool& empty)
{
  uint const pe_size = mesh.topology().comm_size();
  if (pe_size == 1) return;
  if (pe_size & (pe_size - 1))
  {
    propagate_naive(mesh, propagated, empty);
  }
  else
  {
    propagate_hypercube(mesh, propagated, empty);
  }
}
//-----------------------------------------------------------------------------
#ifdef HAVE_MPI
//-----------------------------------------------------------------------------
void DMesh::propagate_naive(Mesh& mesh, Array<Propagation>& propagated, bool& empty)
{
  Comm& comm =  mesh.topology().comm();
  uint pe_rank = mesh.topology().comm_rank();
  uint pe_size = mesh.topology().comm_size();

  // Allocate receive buffer
  int num_prop = propagate.size() * 5;
  int max_prop, recv_count;
  MPI_Allreduce(&num_prop, &max_prop, 1, MPI_INT, MPI_MAX, comm);

  long *recv_buff = new long[max_prop];
  long *send_buff = new long[num_prop];
  long *sp = &send_buff[0];

  for (Array<Propagation>::iterator it = propagate.begin();
       it != propagate.end(); ++it)
  {
    *(sp++) = static_cast<long>(it->first);
    *(sp++) = it->second.mv;
    *(sp++) = it->second.v1;
    *(sp++) = it->second.v2;
    *(sp++) = static_cast<long>(it->second.owner);
  }

  MPI_Status status;
  uint dest, src;

  empty = true;
  for (uint j = 1; j < pe_size; j++)
  {
    src = (pe_rank - j + pe_size) % pe_size;
    dest = (pe_rank + j) % pe_size;

    MPI_Sendrecv(&send_buff[0], num_prop, MPI_LONG, dest, 1, recv_buff,
                 max_prop, MPI_LONG, src, 1, comm, &status);
    MPI_Get_count(&status, MPI_LONG, &recv_count);

    if (recv_count > 0) empty = false;

    dolfin_assert(recv_count % 5 == 0);
    for (int k = 0; k < recv_count; k += 5)
    {

      prop_edge node;
      node.mv = recv_buff[k + 1];
      node.v1 = recv_buff[k + 2];
      node.v2 = recv_buff[k + 3];
      node.owner = static_cast<uint>(recv_buff[k + 4]);

      Propagation prop(static_cast<uint>(recv_buff[k]), node);
      propagated.push_back(prop);
    }

  }

  less_pair comp;
  std::sort(propagated.begin(), propagated.end(), comp);

  short prop, gprop;
  prop = (empty == false);
  MPI_Allreduce(&prop, &gprop, 1, MPI_SHORT, MPI_SUM, comm);
  empty = (gprop == 0);

  delete[] send_buff;
  delete[] recv_buff;
}
//-----------------------------------------------------------------------------
void DMesh::propagate_hypercube(Mesh& mesh, Array<Propagation>& propagated, bool& empty)
{
  Comm& comm =  mesh.topology().comm();
  uint pe_rank = mesh.topology().comm_rank();
  uint pe_size = mesh.topology().comm_size();

  // Allocate receive buffer
  int num_prop = propagate.size() * 5;
  int total_prop, recv_count;
  MPI_Allreduce(&num_prop, &total_prop, 1, MPI_INT, MPI_SUM, comm);

  long *recv_buff = new long[total_prop];
  long *state = new long[total_prop];
  long *sp = &state[0];
  int state_size = 0;

  for (Array<Propagation>::iterator it = propagate.begin();
      it != propagate.end(); ++it)
  {
    *(sp++) = static_cast<long>(it->first);
    *(sp++) = it->second.mv;
    *(sp++) = it->second.v1;
    *(sp++) = it->second.v2;
    *(sp++) = static_cast<long>(it->second.owner);
    state_size += 5;
  }

  MPI_Status status;
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
    dest = pe_rank ^ (D << j);

    MPI_Sendrecv(state, state_size, MPI_LONG, dest, 1, recv_buff, total_prop,
                 MPI_LONG, dest, 1, comm, &status);
    MPI_Get_count(&status, MPI_LONG, &recv_count);

    dolfin_assert(recv_count % 5 == 0);
    for (int k = 0; k < recv_count; k += 5)
    {

      prop_edge node;
      node.mv = recv_buff[k + 1];
      node.v1 = recv_buff[k + 2];
      node.v2 = recv_buff[k + 3];
      node.owner = static_cast<uint>(recv_buff[k + 4]);

      Propagation prop(static_cast<uint>(recv_buff[k]), node);
      propagated.push_back(prop);
    }
    memcpy(sp, recv_buff, recv_count * sizeof(long));
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
void DMesh::renumber_glb(_map<long, uint>& new_global)
{  
  uint const pe_rank = MPI::rank();
  uint const pe_size = MPI::size();
  Array<long> *ghost_buffer = new Array<long>[pe_size];

  uint num_owned = 0;
  uint num_ghost = 0;
  
  for (VertexSet::iterator it = vertices.begin(); it != vertices.end(); ++it)
  {
    DVertex * const dv(*it);
    dolfin_assert(!dv->deleted);

    if (dv->is_shared && dv->owner != pe_rank)
    {
      ghost_buffer[dv->owner].push_back(dv->glb_id);
      num_ghost++;
    }
    else
    {
      num_owned++;
    }
  }

  uint offset = 0;
#if ( MPI_VERSION > 1 )
  MPI_Exscan(&num_owned, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
#else
  MPI_Scan(&num_owned, &offset, 1, MPI_UNSIGNED, MPI_SUM, MPI::DOLFIN_COMM);
  offset -= num_owned;
#endif

  for (VertexSet::iterator it = vertices.begin(); it != vertices.end(); ++it)
  {
    DVertex * const dv(*it);
    dolfin_assert(!dv->deleted);
   
    if (!dv->is_shared || (dv->is_shared && dv->owner == pe_rank))
    {
      new_global[dv->glb_id] = offset++;
    }
  }

  MPI_Status status;
  Array<uint> send_buffer;
  uint src, dest;
  int recv_count;

  /// @todo Reduce comm...
  uint recv_size_gh = 0;
  for (uint i = 0; i < pe_size; i++)
  {
    uint send_size = ghost_buffer[i].size();
    MPI_Reduce(&send_size, &recv_size_gh, 1, 
	       MPI_UNSIGNED, MPI_SUM, i, MPI::DOLFIN_COMM);
  }

  long *recv_ghost = new long[recv_size_gh];
  uint *recv_buff = new uint[num_ghost];

  
  for(uint j=1; j < pe_size; j++){
    src = (pe_rank - j + pe_size) % pe_size;
    dest = (pe_rank + j) % pe_size;
    
    MPI_Sendrecv(&ghost_buffer[dest][0], ghost_buffer[dest].size(),
                 MPI_LONG, dest, 1, recv_ghost, recv_size_gh,
                 MPI_LONG, src, 1, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_LONG, &recv_count);

    for(int k=0; k < recv_count; k++)
      send_buffer.push_back(new_global[recv_ghost[k]]);

    MPI_Sendrecv(&send_buffer[0], send_buffer.size(), MPI_UNSIGNED, src, 2,
                 recv_buff, num_ghost, MPI_UNSIGNED, dest, 2,
                 MPI::DOLFIN_COMM,&status);
    MPI_Get_count(&status,MPI_UNSIGNED,&recv_count);

    for(int j=0; j < recv_count; j++)
      new_global[ghost_buffer[dest][j]] = recv_buff[j];
    send_buffer.clear();
  }

  delete[] recv_ghost;
  delete[] recv_buff;
  for(uint i = 0; i < pe_size; i++)
    ghost_buffer[i].clear();
  delete[] ghost_buffer;

}
//-----------------------------------------------------------------------------
#else
//-----------------------------------------------------------------------------
void DMesh::propagate_naive(Mesh& mesh, Array<Propagation>& propagated, bool& empty)
{
  error("Rivara needs MPI");
}
//-----------------------------------------------------------------------------
void DMesh::propagate_hypercube(Mesh& mesh, Array<Propagation>& propagated, bool& empty)
{
  error("Rivara needs MPI");
}
//-----------------------------------------------------------------------------
void DMesh::renumber_glb(_map<long, uint>& new_global)
{
  error("Rivara needs MPI");
}
//-----------------------------------------------------------------------------
#endif

} /* namespace dolfin */

