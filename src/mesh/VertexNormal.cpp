// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008-2009.
//
// First added:  2007-05-01
// Last changed: 2009-12-30

#include <dolfin/mesh/VertexNormal.h>

#include <dolfin/main/MPI.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>

#include <map>

#define B(row,col,nrow) ((row) + ((nrow)*(col)))

namespace dolfin
{

//-----------------------------------------------------------------------------
VertexNormal::VertexNormal(VertexNormal& other) :
    mesh_(other.mesh_),
    alpha_max_(0.5 * DOLFIN_PI),
    weighting_(none)
{
  *this = other;
}

//-----------------------------------------------------------------------------
VertexNormal::VertexNormal(Mesh& mesh, Type weight) :
    mesh_(mesh),
    alpha_max_(0.5 * DOLFIN_PI),
    weighting_(weight)
{
  uint const nsdim = mesh.topology().dim();

  for (uint d = 0; d < nsdim; ++d)
  {
    basis_.push_back(new MeshFunction<real> [nsdim]);
    for (uint i = 0; i < nsdim; ++i)
    {
      basis_.back()[i].init(mesh, 0);
    }
  }
  vertex_type_.init(mesh, 0);

  ComputeNormal(mesh);
}

//-----------------------------------------------------------------------------
VertexNormal::~VertexNormal()
{
  Clear();
}

//-----------------------------------------------------------------------------
void VertexNormal::Clear()
{
  while (!basis_.empty())
  {
    delete[] basis_.back();
    basis_.pop_back();
  }
}

//-----------------------------------------------------------------------------
VertexNormal& VertexNormal::operator=(VertexNormal& other)
{
  Clear();

  uint const gdim = mesh_.geometry().dim();

  // Initialize data structures
  vertex_type_.init(mesh_, 0);
  for (uint d = 0; d < gdim; ++d)
  {
    basis_.push_back(new MeshFunction<real> [gdim]);
    for (uint i = 0; i < gdim; ++i)
    {
      basis_.back()[i].init(mesh_, 0);
    }

  }

  // Copy data
  for (VertexIterator v(mesh_); !v.end(); ++v)
  {
    vertex_type_.set(*v, other.vertex_type_.get(*v));
    for (uint d = 0; d < gdim; ++d)
    {
      for (uint i = 0; i < gdim; ++i)
      {
        basis_[d][i].set(*v, other.basis()[d][i].get(*v));
      }
    }
  }

  return *this;
}

//-----------------------------------------------------------------------------
void VertexNormal::ComputeSimpleNormal(Mesh& mesh)
{
  message(1, "VertexNormal: Compute normals");
  mesh.renumber();

  uint const nsdim = mesh.geometry().dim();
  uint rank = dolfin::MPI::processNumber();
  uint pe_size = dolfin::MPI::numProcesses();
  Array<real> *send_buff_type_basis = new Array<real> [pe_size];
  Array<uint> *send_buff_index = new Array<uint> [pe_size];
  _map<uint, bool> used_shared;

  for (MeshSharedIterator s(mesh.distdata(), 0); !s.end(); ++s)
  {
    used_shared[s.index()] = false;
  }

  // Create the boundary mesh and shared entities data structures
  BoundaryMesh boundary(mesh, BoundaryMesh::exterior);
  MeshFunction<uint>* cell_map = boundary.data().meshFunction("cell map");
  MeshFunction<uint>* vertex_map = boundary.data().meshFunction("vertex map");

  //
  if (dolfin::MPI::numProcesses() > 1)
  {
    CacheSharedArea(mesh, boundary);
  }

  //-------------------------------------------------------------------------
  real const cosalpha_max = std::cos(alpha_max_);

  if (boundary.numCells())
  {
    // Important !
    // Initialize the connectivities between vertices and boundary cells.
    uint const boundary_nsdim = boundary.topology().dim();
    boundary.init(0, boundary_nsdim);
    // Initialize the connectivities between facets and boundary cells.
    mesh_.init(boundary_nsdim, nsdim);
    for (VertexIterator boundary_vertex(boundary); !boundary_vertex.end();
        ++boundary_vertex)
    {
      uint boundary_id = vertex_map->get(*boundary_vertex);
      uint global_id = mesh.distdata().get_vertex_global(boundary_id);
      bool const vertex_is_shared = mesh.distdata().is_shared(boundary_id, 0);
      bool const vertex_is_ghosted = mesh.distdata().is_ghost(boundary_id, 0);

      //--- Get number of neighbouring facets ---------------------------------
      //message("Get number of neighbouring facets");
      uint NbNeighCells = 0;
      if (dolfin::MPI::numProcesses() > 1 && vertex_is_shared)
      {
        NbNeighCells = num_neigh_cells_[global_id];
      }
      else
      {
        NbNeighCells = boundary_vertex->numEntities(boundary_nsdim);
      }
      // A vertex should not have zero neighbouring facets
      dolfin_assert(NbNeighCells > 0);

      //--- Get facet data ----------------------------------------------------
      //message("Get facet data");
      uint neighbour = 0;   // number of neighbouring facets
      Array<real> normals;  // facet normals
      Array<real> weights;  // facet weights
      GetLocalFacetsData(nsdim, *boundary_vertex, *cell_map, neighbour, normals,
                         weights);

      // Add neighbouring facet contributions from other processes
      if (dolfin::MPI::numProcesses() > 1 && vertex_is_shared)
      {
        Array<real>& shared_weights = shared_facetweights_block_[global_id];
        weights.insert(weights.end(), shared_weights.begin(),
                       shared_weights.end());

        Array<real>& shared_normals = shared_facetnormals_block_[global_id];
        normals.insert(normals.end(), shared_normals.begin(),
                       shared_normals.end());

        neighbour += shared_weights.size();
      } dolfin_assert(neighbour == NbNeighCells); dolfin_assert(normals.size() == nsdim * weights.size());

      //--- Determine vertex type by discriminating surfaces ------------------
      // Add storage for normals to surfaces
      Array<Array<real> > surface_normals;
      Array<real> surface_totalweights;

      // Initialize set of normals to be tested
      std::set<uint> remaining_normals;
      std::set<uint>::iterator it = remaining_normals.begin();
      for (uint i = 0; i < NbNeighCells; ++i)
      {
        it = remaining_normals.insert(it, i);
      }

      uint vertex_type = 0;
      real cosalpha = 0.0;
      uint ref_facet = 0;
      uint cur_facet = 0;
      real wSx = 0.0;
      real nSx[3] =
        { 0.0 };
      while (remaining_normals.size() > 0)
      {
        it = remaining_normals.begin();
        // Initialize new surface normal and weight with reference
        ref_facet = (*it);
        wSx = weights[ref_facet];
        for (uint d = 0; d < nsdim; ++d)
        {
          nSx[d] = wSx * normals[nsdim * ref_facet + d];
        }

        std::set<uint> used_normals;
        used_normals.insert(used_normals.begin(), ref_facet);
        // Loop through remaining normals indexes
        for (++it; it != remaining_normals.end(); ++it)
        {
          cur_facet = (*it);
          dolfin_assert(weights[cur_facet]>0);

          // Compute the scalar product with the reference normal
          cosalpha = 0.0;
#pragma _CRI novector
          for (uint d = 0; d < nsdim; ++d)
          {
            cosalpha += normals[nsdim * ref_facet + d]
                * normals[nsdim * cur_facet + d];
          }
          if (cosalpha > cosalpha_max)
          {
            // Add contribution to surface normal
            for (uint d = 0; d < nsdim; ++d)
            {
              nSx[d] += weights[cur_facet] * normals[nsdim * cur_facet + d];
            }
            wSx += weights[cur_facet];
            // Take into account that the normal is used
            used_normals.insert(used_normals.begin(), cur_facet);
          }
        }

        // Add unit surface normal to the list of surface normals
        //message("Surface normal %i", vertex_type);
        NormalizeVector(nSx);
        surface_normals.push_back(Array<real>());
        surface_normals.back().assign(&nSx[0], &nSx[3]);
        surface_totalweights.push_back(wSx);

        // Next loop we add a new surface and discriminate again against
        // the remaining normals
        for (it = used_normals.begin(); it != used_normals.end(); ++it)
        {
          remaining_normals.erase(*it);
        }
        ++vertex_type;
      } dolfin_assert(surface_totalweights.size() == vertex_type);

      //--- Compute vertex normal ---------------------------------------------
      // The typical algorithm would be:
      // n_k    = sum_{i=1}^k nS_i
      // tau1_k = |n_{k-1}|^2 nS_k - (n_{k-1} . nS_k ) n_{k-1}
      // tau2_k = n_k ^ tau1_k
      // and such that in 2d tau2_k = ez = (0 , 0, 1)
      real basis_vec[3][3] =
        {
          { 1.0, 0.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 1.0 } };

      if (weighting_ == none)
      {
        surface_totalweights.resize(surface_totalweights.size(), 1.0);
      }
      //message("Compute vertex normal");

      for (uint in = 0; in < nsdim; ++in)
      {
        basis_vec[0][in] = 0.0;
        for (uint s = 0; s < vertex_type; ++s)
        {
          basis_vec[0][in] += surface_totalweights[s] * surface_normals[s][in];
        }
      }

      // Normalize the first vector in basis_vec (which is the normal)
      NormalizeVector(basis_vec[0]);

      // Compute unit tangents from unit outward normal
      // Taken from V. John, J. Comp. and Appl. Math. 2002
      switch (nsdim)
        {
        case 1:
          break;
        case 2:
          basis_vec[1][0] = -basis_vec[0][1];
          basis_vec[1][1] = basis_vec[0][0];
          basis_vec[1][2] = 0.0;
          basis_vec[2][0] = 0.0;
          basis_vec[2][1] = 0.0;
          basis_vec[2][2] = 0.0;
          break;
        case 3:
          if (vertex_type == 1)
          {

            if (std::fabs(basis_vec[0][0]) >= 0.5
                || std::fabs(basis_vec[0][1]) >= 0.5)
            {
              // t11 = n2/n
              basis_vec[1][0] = basis_vec[0][1];
              // t12 = -n1/n
              basis_vec[1][1] = -basis_vec[0][0];
              // t13 = 0
              basis_vec[1][2] = 0.0;

              NormalizeVector(basis_vec[1]);

              // t21 = -t12*n3
              basis_vec[2][0] = -basis_vec[1][1] * basis_vec[1][2];
              // t22 = t11*n3
              basis_vec[2][1] = basis_vec[1][0] * basis_vec[1][2];
              // t23 = t12*n1 - t11*n2
              basis_vec[2][2] = basis_vec[1][1] * basis_vec[0][0]
                  - basis_vec[1][0] * basis_vec[0][1];

              // basis_vec[2] is de facto normalized
            }
            else
            {
              // t11 = 0
              basis_vec[1][0] = 0.0;
              // t12 = -n3/n
              basis_vec[1][1] = -basis_vec[0][2];
              // t13 = n2/n
              basis_vec[1][2] = basis_vec[0][1];

              NormalizeVector(basis_vec[1]);

              // t21 = t13*n2 - t12*n3
              basis_vec[2][0] = basis_vec[1][2] * basis_vec[0][1]
                  - basis_vec[1][1] * basis_vec[0][2];
              // t22 = -t13*n1
              basis_vec[2][1] = -basis_vec[1][2] * basis_vec[0][0];
              // t23 = t12*n1
              basis_vec[2][2] = basis_vec[1][1] * basis_vec[0][0];

              // basis_vec[2] is de facto normalized
            }
          }
          else
          {
            uint const Sk = vertex_type - 1;
            basis_vec[2][0] = basis_vec[0][1] * surface_normals[Sk][2]
                - surface_normals[Sk][1] * basis_vec[0][2];
            basis_vec[2][1] = basis_vec[0][2] * surface_normals[Sk][0]
                - surface_normals[Sk][2] * basis_vec[0][0];
            basis_vec[2][2] = basis_vec[0][0] * surface_normals[Sk][1]
                - surface_normals[Sk][0] * basis_vec[0][1];

            NormalizeVector(basis_vec[2]);

            basis_vec[1][0] = basis_vec[2][1] * basis_vec[0][2]
                - basis_vec[2][2] * basis_vec[0][1];
            basis_vec[1][1] = basis_vec[2][2] * basis_vec[0][0]
                - basis_vec[2][0] * basis_vec[0][2];
            basis_vec[1][2] = basis_vec[2][0] * basis_vec[0][1]
                - basis_vec[2][1] * basis_vec[0][0];

            // basis_vec[2] is de facto normalized
          }
          break;
        default:
          break;
        }

      //--- Distribute vertex normals -----------------------------------------
      // Prepare data structures
      if (vertex_is_ghosted)
      {
        uint owner = mesh.distdata().get_owner(boundary_id, 0);
        send_buff_type_basis[owner].push_back((double) vertex_type);
        for (uint basisidx = 0; basisidx < nsdim * nsdim; ++basisidx)
        {
          for (uint in = 0; in < nsdim; ++in)
          {
            send_buff_type_basis[owner].push_back(basis_vec[basisidx][in]);
          }
        }
        send_buff_index[owner].push_back(global_id);
      }
      else
      {
        vertex_type_.set(boundary_id, vertex_type);
        if (vertex_type == 0)
        {
          error("Surface multiplicity is equal to zero");
        }
        for (uint basisvec = 0; basisvec < nsdim; ++basisvec)
        {
          for (uint in = 0; in < nsdim; ++in)
          {
            basis_[basisvec][in].set(boundary_id, basis_vec[basisvec][in]);
          }
        }
        used_shared[boundary_id] = true;
      }
    }
  }

  // Synchronize basis vectors, vertex_types across processes
  if (dolfin::MPI::numProcesses() > 1)
  {
#ifdef HAVE_MPI
    MPI_Status status;
    uint src = 0;
    uint dest = 0;
    int recv_count = 0;
    int recv_size = 0;
    int send_size = 0;
    int recv_count_data = 0;

    // Collect data size
    for (uint i = 0; i < pe_size; i++)
    {
      send_size = send_buff_index[i].size();
      MPI_Reduce(&send_size, &recv_count, 1, MPI_INT, MPI_SUM, i,
                 dolfin::MPI::DOLFIN_COMM);

      send_size = send_buff_type_basis[i].size();
      MPI_Reduce(&send_size, &recv_count_data, 1, MPI_INT, MPI_SUM, i,
                 dolfin::MPI::DOLFIN_COMM);
    }

    // Storage is (vertex_type (size = 1), basis (size = nsdim*nsdim))
    uint data_alignment = 1 + nsdim * nsdim;
    uint *recv_index = new uint[recv_count];
    real *recv_type = new real[recv_count_data];

    for (uint i = 1; i < pe_size; i++)
    {
      src = (rank - i + pe_size) % pe_size;
      dest = (rank + i) % pe_size;

      MPI_Sendrecv(&send_buff_index[dest][0], send_buff_index[dest].size(),
                   MPI_UNSIGNED, dest, 0, recv_index, recv_count, MPI_UNSIGNED,
                   src, 0, dolfin::MPI::DOLFIN_COMM, &status);

      MPI_Sendrecv(&send_buff_type_basis[dest][0],
                   send_buff_type_basis[dest].size(), MPI_DOUBLE, dest, 1,
                   recv_type, recv_count_data, MPI_DOUBLE, src, 1,
                   dolfin::MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &recv_size);
      // Insert check if value assigned
      uint idx = 0;
      // Data alignment is n_tau + 1
      for (int j = 0; j < recv_size; j += data_alignment, ++idx)
      {
        uint index = mesh.distdata().get_vertex_local(recv_index[idx]);
        if (!used_shared[index])
        {
          vertex_type_.set(index, (uint) recv_type[j]);

          uint offset = j + 1;
          for (uint basisvec = 0; basisvec < nsdim; ++basisvec)
          {
            for (uint in = 0; in < nsdim; ++in)
            {
              basis_[basisvec][in].set(
                  index, recv_type[offset + basisvec * nsdim + in]);
            }
          }
          used_shared[index] = true;
        }
      }
    }

    delete[] recv_index;
    delete[] recv_type;
#endif
  }

//--- Cleanup
  for (uint i = 0; i < pe_size; ++i)
  {
    send_buff_type_basis[i].clear();
    send_buff_index[i].clear();
  }
  delete[] send_buff_type_basis;
  delete[] send_buff_index;

  num_neigh_cells_.clear();
  shared_offsetidx_.clear();

  for (std::map<uint, Array<real> >::iterator it =
      shared_facetnormals_block_.begin();
      it != shared_facetnormals_block_.end(); ++it)
  {
    it->second.clear();
  }
  for (std::map<uint, Array<real> >::iterator it =
      shared_facetweights_block_.begin();
      it != shared_facetweights_block_.end(); ++it)
  {
    it->second.clear();
  }
  shared_facetnormals_block_.clear();
  shared_facetweights_block_.clear();
}

//-----------------------------------------------------------------------------
void VertexNormal::GetLocalFacetsData(uint const& gdim, Vertex& vertex,
                                      MeshFunction<uint>& cell_map,
                                      uint& nb_neigh, Array<real>& normals,
                                      Array<real>& weights)
{
  //
  for (CellIterator boundary_cell(vertex); !boundary_cell.end();
      ++boundary_cell)
  {
    // Create mesh facet corresponding to boundary cell
    Facet mesh_facet(mesh_, cell_map.get(*boundary_cell));
    dolfin_assert(mesh_facet.numEntities(gdim) == 1);

    // Get cell to which facet belongs (pick first, there is only one)
    Cell mesh_cell(mesh_, mesh_facet.entities(gdim)[0]);

    // Get local index of facet with respect to the cell
    uint local_facet = mesh_cell.index(mesh_facet);

    // Add facet weight
    switch (weighting_)
      {
      case none:
        weights.push_back(1.0);
        break;
      case facet:
        // Compute the measure of the facet
        weights.push_back(boundary_cell->volume());
        break;
      case cell:
        // Compute the measure of the cell
        weights.push_back(mesh_cell.volume());
        break;
      }

    // Add facet normal components
    for (uint d = 0; d < gdim; ++d)
    {
      normals.push_back(mesh_cell.normal(local_facet, d));
    }

    ++nb_neigh;
  }
}

//-----------------------------------------------------------------------------
void VertexNormal::CacheSharedArea(Mesh& mesh, BoundaryMesh& boundary)
{
#ifdef HAVE_MPI
  uint const nsdim = mesh.geometry().dim();

  int rank = dolfin::MPI::processNumber();
  int pe_size = dolfin::MPI::numProcesses();

  // Send buff for global indices of shared vertices
  Array<uint> sendbuff_global_vert_indices;
  // Send buff for facets normals associated with shared vertices
  Array<real> sendbuff_facetnormals;
  // Send buff for weights of shared vertices
  Array<real> sendbuff_facetweights;

  // Send buff for offset indices packed by triplet for each shared vertex
  // (NbNeighbouringCells, FacetNormalsOffset, FacetWeightsOffset )
  Array<uint> sendbuff_offset_indices;
  vertex_offset_ = 0;
  facetnormals_offset_ = 0;
  facetweights_offset_ = 0;

  // Map global id of boundary vertices to boundary marker
  std::map<uint, bool> GlobalIdOnBoundary;

  uint SharedVertexCount = 0;
  uint SharedMeshFacetCount = 0;

  // Computation of normals to the boundary vertices shared between processes
  if (boundary.numCells())
  {
    uint NbNeighCells = 0;
    MeshFunction<uint> *cell_map = boundary.data().meshFunction("cell map");
    MeshFunction<uint> *vertex_map = boundary.data().meshFunction("vertex map");
    for (VertexIterator boundary_vertex(boundary); !boundary_vertex.end();
        ++boundary_vertex, SharedMeshFacetCount += NbNeighCells)
    {
      uint boundary_id = vertex_map->get(*boundary_vertex);
      uint global_id = mesh.distdata().get_vertex_global(boundary_id);

      // Mark vertex as part of the boundary
      GlobalIdOnBoundary[global_id] = true;

      // Cache number of neighbouring cells for each shared vertex
      NbNeighCells = 0;

      bool const vertex_is_shared = mesh.distdata().is_shared(boundary_id, 0);

      if (vertex_is_shared)
      {
        ++SharedVertexCount;
        sendbuff_global_vert_indices.push_back(global_id);

        // Get local facet data
        GetLocalFacetsData(nsdim, *boundary_vertex, *cell_map, NbNeighCells,
                           sendbuff_facetnormals, sendbuff_facetweights);

        num_neigh_cells_[global_id] = NbNeighCells;

        sendbuff_offset_indices.push_back(NbNeighCells);
        sendbuff_offset_indices.push_back(facetnormals_offset_);
        sendbuff_offset_indices.push_back(facetweights_offset_);

        // Init datastructures for shared data
        shared_offsetidx_[global_id] = vertex_offset_;

        // Update shared vertex offset
        vertex_offset_ += nsdim;
        // Padding for NbNeighCells FacetNormals
        facetnormals_offset_ += nsdim * NbNeighCells;
        // Padding for NbNeighCells FacetWeights
        facetweights_offset_ += NbNeighCells;
      }
    }
  }

// Exchange values
  MPI_Status status;
  uint src;
  uint dest;

  int sh_vertidx_count = sendbuff_global_vert_indices.size();
  int sh_facetnormals_count = sendbuff_facetnormals.size();
  int sh_facetweights_count = sendbuff_facetweights.size();

  dolfin_assert(sh_vertidx_count == (int) SharedVertexCount);
  dolfin_assert(sh_facetnormals_count == (int) SharedMeshFacetCount* (int) nsdim);
  dolfin_assert(sh_facetweights_count == (int) SharedMeshFacetCount);

  int recv_size_vertidx;
  int recv_size_facetnormals;
  int recv_size_facetweights;

  int recv_count;

  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
  MPI_Allreduce(&sh_vertidx_count, &recv_size_vertidx, 1, MPI_INT, MPI_MAX,
                dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&sh_facetweights_count, &recv_size_facetweights, 1, MPI_INT,
                MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&sh_facetnormals_count, &recv_size_facetnormals, 1, MPI_INT,
                MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  uint *recv_vertidx = new uint[recv_size_vertidx];
  uint *recv_offsetidx = new uint[offsetidx_padding_ * recv_size_vertidx];

  real *recv_facetnormals = new real[recv_size_facetnormals];
  real *recv_facetweights = new real[recv_size_facetweights];

// For each process
  for (int proc = 1; proc < pe_size; ++proc)
  {
    src = (rank - proc + pe_size) % pe_size;
    dest = (rank + proc) % pe_size;

    MPI_Sendrecv(&sendbuff_global_vert_indices[0], sh_vertidx_count,
                 MPI_UNSIGNED, src, 1, recv_vertidx, recv_size_vertidx,
                 MPI_UNSIGNED, dest, 1, dolfin::MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    MPI_Sendrecv(&sendbuff_facetnormals[0], sh_facetnormals_count, MPI_DOUBLE,
                 src, 1, recv_facetnormals, recv_size_facetnormals, MPI_DOUBLE,
                 dest, 1, dolfin::MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&sendbuff_facetweights[0], sh_facetweights_count, MPI_DOUBLE,
                 src, 1, recv_facetweights, recv_size_facetweights, MPI_DOUBLE,
                 dest, 1, dolfin::MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&sendbuff_offset_indices[0],
                 (offsetidx_padding_ * sh_vertidx_count), MPI_UNSIGNED, src, 1,
                 recv_offsetidx, offsetidx_padding_ * recv_size_vertidx,
                 MPI_UNSIGNED, dest, 1, dolfin::MPI::DOLFIN_COMM, &status);

    // Index for vertex
    uint vertidx = 0;
    // Index for boundary cells == global mesh facets
    uint facetidx = 0;
    // Number of neighbouring boundary cells
    uint nbneighcells = 0;
    // Offsets
    uint facetnoffset = 0;
    uint weightoffset = 0;
    for (int i = 0; i < recv_count; i++)
    {
      uint glb_index = recv_vertidx[i];
      if (mesh.distdata().has_global(glb_index, 0)
          && GlobalIdOnBoundary.count(glb_index) > 0
          && GlobalIdOnBoundary[glb_index])
      {
        // Get alignment and offsets
        nbneighcells = recv_offsetidx[facetidx];
        facetnoffset = recv_offsetidx[facetidx + 1];
        weightoffset = recv_offsetidx[facetidx + 2];

        // Update number of neighbouring boundary cells
        num_neigh_cells_[glb_index] += nbneighcells;

        // Add corresponding facet normals
        for (uint k = 0; k < nsdim * nbneighcells; ++k)
        {
          shared_facetnormals_block_[glb_index].push_back(
              recv_facetnormals[facetnoffset++]);
        }

        // Add corresponding facet weights
        for (uint k = 0; k < nbneighcells; ++k)
        {
          shared_facetweights_block_[glb_index].push_back(
              recv_facetweights[weightoffset++]);
        }
      }

      // Increase data offsets
      vertidx += nsdim;
      facetidx += offsetidx_padding_;
    }

  }
  delete[] recv_vertidx;
  delete[] recv_facetnormals;
  delete[] recv_offsetidx;
  delete[] recv_facetweights;
#endif
}

//-----------------------------------------------------------------------------
void VertexNormal::ComputeNormal(Mesh& mesh)
{
  mesh.renumber();
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  Array<real> *send_buff_type = new Array<real> [pe_size];
  Array<uint> *send_buff_index = new Array<uint> [pe_size];
  _map<uint, bool> used_shared;

  for (MeshSharedIterator s(mesh.distdata(), 0); !s.end(); ++s)
  {
    used_shared[s.index()] = false;
  }

  int nsdim = mesh.topology().dim();

  // Iterate over all cells in the boundary mesh
  BoundaryMesh boundary(mesh, BoundaryMesh::exterior);
  MeshFunction<uint>* cell_map = boundary.data().meshFunction("cell map");
  MeshFunction<uint>* vertex_map = boundary.data().meshFunction("vertex map");
  if (boundary.numCells() > 0)
  {
    Progress p_boundary("Assembling boundary contributions",
                        boundary.numCells());
  }

  real normal_vec[3];
  real normal_vec1[3];
  real tau_vec[2];
  real tau_vec_1[3];
  real tau_vec_2[3];

  real *n_block = 0;
  real *ns1_block = 0;
  real *ns2_block = 0;
  real *area_block = 0;

  if (MPI::numProcesses() > 1) CacheSharedArea(mesh, boundary);

  // Computation of normals to the boundary vertices
  if (boundary.numCells() > 0)
  {
    for (VertexIterator n(boundary); !n.end(); ++n)
    {

      int id = vertex_map->get(*n);

      int Nrow = 0;
      for (CellIterator boundary_cell(*n); !boundary_cell.end();
          ++boundary_cell)
        Nrow++;

      // Add storage for shared vertices cell normals
      if (MPI::numProcesses() > 1) if (mesh.distdata().is_shared(id, 0)) Nrow =
          num_neigh_cells_[mesh.distdata().get_global(id, 0)];

      int Ncol = 3;

      n_block = new real[Nrow * Ncol]; // contains of all normals from the surfaces

      int irow = 0;
      area_block = new real[Nrow];  // contains of all areas from the elements

      real sum_area = 0.0;         // sum of all areas of the elements

      // in a case of node on the edge
      real sum_area_s1 = 0.0; // sum of all areas of the elements which belong to S1
      real sum_area_s2 = 0.0; // sum of all areas of the elements which belong to S2

      for (int i = 0; i < nsdim; i++)
      {
        normal_vec[i] = 0.0;
        normal_vec1[i] = 0.0;
      }

      //for (CellIterator boundary_cell(boundary); !boundary_cell.end(); ++boundary_cell)
      for (CellIterator boundary_cell(*n); !boundary_cell.end();
          ++boundary_cell)
      {
        // Create mesh facet corresponding to boundary cell
        Facet mesh_facet(mesh, cell_map->get(*boundary_cell));
        dolfin_assert(mesh_facet.numEntities(nsdim) == 1);

        // Get cell to which facet belongs (pick first, there is only one)
        Cell mesh_cell(mesh, mesh_facet.entities(mesh.topology().dim())[0]);
        // Get local index of facet with respect to the cell
        uint local_facet = mesh_cell.index(mesh_facet);

        // Compute the lenght/area of the cell
        real area = mesh_cell.volume();
        area_block[irow] = area;
        // Get sum of the length/area of all neighbouring cells
        sum_area += area;

        // Compute a normal to the boundary
        normal_vec[0] += area * mesh_cell.normal(local_facet, 0);
        normal_vec[1] += area * mesh_cell.normal(local_facet, 1);
        if (nsdim == 3) normal_vec[2] += area
            * mesh_cell.normal(local_facet, 2);

        for (int j = 0; j < Ncol; j++)
        {
          n_block[B(irow, j, Nrow)] = mesh_cell.normal(local_facet, j);
        }
        irow++;
      }

      if (MPI::numProcesses() > 1) if (mesh.distdata().is_shared(id, 0))
      {
        uint glb_index = mesh.distdata().get_global(id, 0);
        //    sum_area += shared_area[glb_index];
        normal_vec[0] += shared_normal[shared_offsetidx_[glb_index]];
        normal_vec[1] += shared_normal[shared_offsetidx_[glb_index] + 1];
        if (nsdim == 3) normal_vec[2] += shared_normal[shared_offsetidx_[glb_index]
            + 2];

        Array<real> tmp = shared_facetweights_block_[glb_index];
        std::vector<real>::iterator iter;

        int iirow = irow;
        for (iter = tmp.begin(); iter != tmp.end(); iter++)
        {
          area_block[iirow] = *iter;
          sum_area += area_block[iirow];
          iirow++;
        }

        tmp = shared_facetnormals_block_[glb_index];
        for (iter = tmp.begin(); iter != tmp.end(); iter += 3)
        {
          for (int j = 0; j < Ncol; j++)
            n_block[B(irow,j,Nrow)] = *(iter + j);
          irow++;
        }
      }

      ns1_block = new real[Nrow * Ncol]; // contains of all normals from the case S1
      ns2_block = new real[Nrow * Ncol]; // contains of all normals from the case S2

      real n1[3];
      int n_type = 1;                 // n_type is a node or vertex type:
      // n_type = 1 vector in S1,
      // n_type = 2 in S2,
      // n_type = 3 in S3
      // S1 - vertex lies on the surface,
      // S2 - edge, S3 - corner
      real alpha = 0.0;        // needed for the angle between the vectors

      // take one normal from the matrix of normals n1_block
      for (int j = 0; j < Ncol; j++)
        n1[j] = n_block[B(0,j,Nrow)];

      real alpha_max = DOLFIN_PI / 6;
      int Nrow_ns1 = 0;  // it will be number of rows for the ns1_block
      int Nrow_ns2 = 0;  // it will be number of rows for the ns2_block

      // find all vector which are from S2 - edge

      for (int i = 0; i < Nrow; i++)
      {
        real csum1 = 0.0;
        real csum2 = 0.0;

        for (int in = 0; in < nsdim; in++)
        {
          csum1 += sqr(n1[in]);
          csum2 += sqr(n_block[B(i, in, Nrow)]);
        }

        real n_norm = sqrt(csum1);
        real nj_norm = sqrt(csum2);

        real csum3 = 0.0;
        for (int in = 0; in < nsdim; in++)
          csum3 += n1[in] * n_block[B(i,in,Nrow)];

        alpha = acos(csum3 / (n_norm * nj_norm));

        if (nsdim == 2 && fabs(alpha) > alpha_max) n_type = 2;

        if (nsdim == 3)
        {

          if (fabs(alpha) > alpha_max)
          {
            n_type = 2;

            // save normals to the block ns2_block
            for (int j = 0; j < Ncol; j++)
            {
              ns2_block[B(Nrow_ns2,j,Nrow)] = n_block[B(i,j,Nrow)];
            }

            // calculate the area of all elements belongs ot S2
            sum_area_s2 += area_block[i];
            Nrow_ns2++;
          }
          else
          {
            // save normals to the block ns1_block
            for (int j = 0; j < Ncol; j++)
              ns1_block[B(Nrow_ns1, j,Nrow)] = n_block[B(i,j,Nrow)];

            // calculate the area of all elements belongs ot S2
            sum_area_s1 += area_block[i];
            Nrow_ns1++;
          }
        }
      }  // end of for

      if (nsdim == 3)
      {
        // go through all normals from S2 and find normals from S3 if any
        if (n_type == 2)
        {
          for (int j = 0; j < Ncol; j++)
            n1[j] = ns2_block[B(0,j,Nrow)];

          // find all vector which are from S3 - corner
          for (int i = 0; i < Nrow_ns2; i++)
          {
            real n_norm = sqrt(sqr(n1[0]) + sqr(n1[1]) + sqr(n1[2]));
            real nj_norm = sqrt(
                sqr(ns2_block[B(i,0,Nrow)]) + sqr(ns2_block[B(i,1,Nrow)])
                    + sqr(ns2_block[B(i,2,Nrow)]));
            alpha = acos(
                (n1[0] * ns2_block[B(i,0,Nrow)] + n1[1] * ns2_block[B(i,1,Nrow)]
                    + n1[2] * ns2_block[B(i,2,Nrow)]) / n_norm * nj_norm);
            if (fabs(alpha) > alpha_max)
            {
              n_type = 3;
            }
          }
        }

        // Compute a normal to the vertex
        if (n_type == 1)
        {  // the vertex lies on the surface, case S1
          for (int in = 0; in < nsdim; in++)
          {
            normal_vec[in] /= sum_area;
          }
        }
        else if (n_type == 2)
        {  // the vertex lies on the edge, case S2
           // now we use the matrix of normals for S2: ns1_block and ns2_block
           // note ns1_block is used for the first normal and
           // ns2_block is for the second normal at the same edge
           // here we take the second normal instead of one tangental vector

          // compute the first normal from the normals of S1
          for (int in = 0; in < nsdim; in++)
            normal_vec[in] = 0.0;

          for (int i = 0; i < Nrow_ns1; i++)
            for (int in = 0; in < nsdim; in++)
              normal_vec[in] += ns1_block[B(i,in,Nrow)];

          for (int in = 0; in < nsdim; in++)
            normal_vec[in] /= sum_area_s1;

          // compute the second normal from the normals of S2
          for (int i = 0; i < Nrow_ns2; i++)
            for (int in = 0; in < nsdim; in++)
              normal_vec1[in] += ns2_block[B(i,in,Nrow)];

          for (int in = 0; in < nsdim; in++)
            normal_vec1[in] /= sum_area_s2;
        }
      }

      // Find the length of the normal vector and normal_vecize it
      real l = 0.0;
      for (int in = 0; in < nsdim; in++)
        l += sqr(normal_vec[in]);
      l = sqrt(l);

      for (int in = 0; in < nsdim; in++)
      {
        normal_vec[in] /= l;

        if (fabs(normal_vec[in]) < DOLFIN_EPS) // check if the component is too small replace it by zero
        normal_vec[in] = 0.0;
      }

      if (nsdim == 3 && n_type == 2)
      {  // normalize normal1 if the case S2 exist
        l = 0.0;
        for (int in = 0; in < nsdim; in++)
          l += sqr(normal_vec1[in]);
        l = sqrt(l);

        for (int in = 0; in < nsdim; in++)
        {
          normal_vec1[in] /= l;

          if (fabs(normal_vec1[in]) < DOLFIN_EPS) // check if the component is too small replace it by zero
          normal_vec1[in] = 0.0;
        }

      }

      if (!mesh.distdata().is_ghost(id, 0))
      {
        for (int in = 0; in < nsdim; in++)
          basis_[0][in].set(id, normal_vec[in]);
        vertex_type_.set(id, n_type);
        used_shared[id] = true;
      }
      else
      {
        uint owner = mesh.distdata().get_owner(id, 0);
        send_buff_type[owner].push_back(n_type);
        for (int in = 0; in < nsdim; in++)
          send_buff_type[owner].push_back(normal_vec[in]);
        send_buff_index[owner].push_back(mesh.distdata().get_global(id, 0));
      }

      if (nsdim == 2)
      {
        // 2D case
        tau_vec[0] = normal_vec[1];
        tau_vec[1] = -normal_vec[0];

        // Find the length of the tau_vec vector and normalize it
        l = 0.0;
        for (int in = 0; in < nsdim; in++)
          l += sqr(tau_vec[in]);
        l = sqrt(l);

        for (int in = 0; in < nsdim; in++)
        {
          tau_vec[in] /= l;

          // check if the component is too small replace it by zero
          if (fabs(tau_vec[in]) < DOLFIN_EPS) tau_vec[in] = 0.0;
        }
        if (!mesh.distdata().is_ghost(id, 0))
        {
          for (int in = 0; in < nsdim; in++)
            basis_[1][in].set(id, tau_vec[in]);
        }
        else
        {
          uint owner = mesh.distdata().get_owner(id, 0);
          send_buff_type[owner].push_back(tau_vec[0]);
          send_buff_type[owner].push_back(tau_vec[1]);
        }
      }
      else
      {
        // 3D case
        real norm_inv = 0.0;
        if (fabs(normal_vec[0]) >= 0.5 || fabs(normal_vec[1]) >= 0.5)
        {
          norm_inv = 1 / sqrt(sqr(normal_vec[0] + sqr(normal_vec[1])));
          if (n_type == 2)
          {
            tau_vec_1[0] = normal_vec1[0];
            tau_vec_1[1] = normal_vec1[1];
            tau_vec_1[2] = normal_vec1[2];

            // cross product of tau_1 and normal_vec
            tau_vec_2[0] = normal_vec[1] * tau_vec_1[2]
                - tau_vec_1[1] * normal_vec[2];
            tau_vec_2[1] = normal_vec[2] * tau_vec_1[0]
                - tau_vec_1[2] * normal_vec[0];
            tau_vec_2[2] = normal_vec[0] * tau_vec_1[1]
                - tau_vec_1[0] * normal_vec[1];
          }
          else
          {
            tau_vec_1[0] = normal_vec[1] * norm_inv;
            tau_vec_1[1] = -normal_vec[0] * norm_inv;
            tau_vec_1[2] = 0.0;

            tau_vec_2[0] = -tau_vec_1[1] * normal_vec[2];
            tau_vec_2[1] = tau_vec_1[0] * normal_vec[2];
            tau_vec_2[2] = tau_vec_1[1] * normal_vec[0]
                - tau_vec_1[0] * normal_vec[1];
          }
          norm_inv = 0.0;
        }
        else
        {
          norm_inv = 1 / sqrt(sqr(normal_vec[1] + sqr(normal_vec[2])));
          if (n_type == 2)
          {
            tau_vec_1[0] = normal_vec1[0];
            tau_vec_1[1] = normal_vec1[1];
            tau_vec_1[2] = normal_vec1[2];

            // cross product of tau_vec_1 and normal_vec
            tau_vec_2[0] = normal_vec[1] * tau_vec_1[2]
                - tau_vec_1[1] * normal_vec[2];
            tau_vec_2[1] = normal_vec[2] * tau_vec_1[0]
                - tau_vec_1[2] * normal_vec[0];
            tau_vec_2[2] = normal_vec[0] * tau_vec_1[1]
                - tau_vec_1[0] * normal_vec[1];
          }
          else
          {
            tau_vec_1[0] = 0.0;
            tau_vec_1[1] = -normal_vec[2] * norm_inv;
            tau_vec_1[2] = normal_vec[1] * norm_inv;

            tau_vec_2[0] = tau_vec_1[2] * normal_vec[1]
                - tau_vec_1[1] * normal_vec[2];
            tau_vec_2[1] = -tau_vec_1[2] * normal_vec[0];
            tau_vec_2[2] = tau_vec_1[1] * normal_vec[0];
          }
          norm_inv = 0.0;
        }

        // Find the length of the tau_vec_1 vector and normalize it
        l = 0.0;
        for (int in = 0; in < nsdim; in++)
          l += sqr(tau_vec_1[in]);
        l = sqrt(l);
        for (int in = 0; in < nsdim; in++)
        {
          tau_vec_1[in] /= l;

          // check if the component is too small replace it by zero
          if (fabs(tau_vec_1[in]) < DOLFIN_EPS) tau_vec_1[in] = 0.0;
        }

        // Find the length of the tau_vec_2 vector and normalize it
        l = 0.0;
        for (int in = 0; in < nsdim; in++)
        {
          l += sqr(tau_vec_2[in]);
        }
        l = sqrt(l);
        for (int in = 0; in < nsdim; in++)
        {
          tau_vec_2[in] /= l;

          if (fabs(tau_vec_2[in]) < DOLFIN_EPS) // check if the component is too small replace it by zero
          tau_vec_2[in] = 0.0;
        }

        if (!mesh.distdata().is_ghost(id, 0))
        {

          for (int in = 0; in < nsdim; in++)
          {
            basis_[1][in].set(id, tau_vec_1[in]);
            basis_[2][in].set(id, tau_vec_2[in]);
          }
        }
        else
        {
          uint owner = mesh.distdata().get_owner(id, 0);
          send_buff_type[owner].push_back(tau_vec_1[0]);
          send_buff_type[owner].push_back(tau_vec_1[1]);
          send_buff_type[owner].push_back(tau_vec_1[2]);

          send_buff_type[owner].push_back(tau_vec_2[0]);
          send_buff_type[owner].push_back(tau_vec_2[1]);
          send_buff_type[owner].push_back(tau_vec_2[2]);
        }

      }

      delete[] area_block;
      delete[] n_block;
      delete[] ns1_block;
      delete[] ns2_block;

    }
  }

#ifdef HAVE_MPI
  if (MPI::numProcesses() > 1)
  {
    MPI_Status status;
    uint src, dest;
    int recv_count, recv_size, send_size, recv_count_data;
    recv_count_data = 0;  // == 0
    for (uint i = 0; i < pe_size; i++)
    {
      send_size = send_buff_index[i].size();
      MPI_Reduce(&send_size, &recv_count, 1, MPI_INT, MPI_SUM, i,
                 dolfin::MPI::DOLFIN_COMM);

      send_size = send_buff_type[i].size();
      MPI_Reduce(&send_size, &recv_count_data, 1, MPI_INT, MPI_SUM, i,
                 dolfin::MPI::DOLFIN_COMM);
    }

    uint n_tau = nsdim + nsdim;
    if (nsdim == 3) n_tau = nsdim * 2 + nsdim;
    uint *recv_index = new uint[recv_count];
    real *recv_type = new real[recv_count_data];

    for (uint i = 1; i < pe_size; i++)
    {
      src = (rank - i + pe_size) % pe_size;
      dest = (rank + i) % pe_size;

      MPI_Sendrecv(&send_buff_index[dest][0], send_buff_index[dest].size(),
                   MPI_UNSIGNED, dest, 0, recv_index, recv_count, MPI_UNSIGNED,
                   src, 0, dolfin::MPI::DOLFIN_COMM, &status);

      MPI_Sendrecv(&send_buff_type[dest][0], send_buff_type[dest].size(),
                   MPI_DOUBLE, dest, 1, recv_type, recv_count_data, MPI_DOUBLE,
                   src, 1, dolfin::MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &recv_size);
      // Insert check if value assigned
      uint idx = 0;
      for (int j = 0; j < recv_size; j += (n_tau + 1))
      {
        uint index = mesh.distdata().get_local(recv_index[idx], 0);
        if (!used_shared[index])
        {
          vertex_type_.set(index, (int) recv_type[j]);
          //    node_vec_arr[l2vidx[index]] =  recv_type[j];
          uint offset = 0;
          if (nsdim == 2)
          {
            for (int in = 0; in < nsdim; in++)
            {
              offset = j + in + 1;
              basis_[0][in].set(index, recv_type[offset]);
              basis_[1][in].set(index, recv_type[offset + nsdim]);
            }
          }
          else
          {
            for (int in = 0; in < nsdim; in++)
            {
              offset = j + in + 1;
              basis_[0][in].set(index, recv_type[offset]);
              basis_[1][in].set(index, recv_type[offset + nsdim]);
              basis_[2][in].set(index, recv_type[offset + 2 * nsdim]);
            }
          }
          used_shared[index] = true;
        }
        idx++;
      }
    }

    delete[] recv_index;
    delete[] recv_type;
  }

  for (uint i = 0; i < pe_size; i++)
  {
    send_buff_type[i].clear();
    send_buff_index[i].clear();
  }
  delete[] send_buff_type;
  delete[] send_buff_index;

  shared_offsetidx_.clear();
  num_neigh_cells_.clear();
  shared_normal.clear();

  for (std::map<uint, Array<real> >::iterator it = shared_facetnormals_block_.begin();
      it != shared_facetnormals_block_.end(); it++)
    it->second.clear();
  for (std::map<uint, Array<real> >::iterator it = shared_facetweights_block_.begin();
      it != shared_facetweights_block_.end(); it++)
    it->second.clear();

  shared_facetnormals_block_.clear();
  shared_facetweights_block_.clear();
#endif // HAVE_MPI
}
//-----------------------------------------------------------------------------

}

