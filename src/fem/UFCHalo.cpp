// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/fem/UFCHalo.h>

#include <dolfin/common/AdjacentMapping.h>
#include <dolfin/common/Array.h>
#include <dolfin/fem/Coefficient.h>
#include <dolfin/fem/DofMapSet.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Facet.h>

#include <string>

namespace dolfin
{

//-----------------------------------------------------------------------------
UFCHalo::UFCHalo(UFC& ufc, Array<Coefficient*> const& coefficients,
                 DofMapSet const& dof_map_set) :
    cell0(ufc.cell0),
    cell1(ufc.cell1),
    macro_w(ufc.macro_w),
    facet0(ufc.facet0),
    facet1(ufc.facet1),
    macro_dofs(ufc.macro_dofs),
    ufc_(ufc),
    mesh_(*const_cast<Mesh *>(ufc.mesh.mesh)),
    coefficients_(coefficients),
    dof_map_set_(dof_map_set),
    rank_offsets_(),
    facet_map_(),
    r_packet_size_(0),
    r_data0_(nullptr),
    r_data1_(nullptr),
    u_packet_size_(0),
    u_data0_(nullptr),
    u_data1_(nullptr)
{
  //
  init();
}

//-----------------------------------------------------------------------------
UFCHalo::~UFCHalo()
{
  clear();
}

//-----------------------------------------------------------------------------
void UFCHalo::init()
{
  // Early exit as nothing has to be done.
  if (!mesh_.is_distributed())
  {
    return;
  }

  // Clear data structures and define data padding
  clear();
  uint const gdim = mesh_.geometry_dimension();
  uint const facet_dim = mesh_.type().facet_dim();
  DistributedData& distdata = mesh_.distdata()[facet_dim];

  //
  uint const num_cell_vertices = mesh_.type().num_entities(0);
  uint const coordinates_data_size = num_cell_vertices * gdim;
  uint coefficients_data_size = 0;
  for (uint i = 0; i < ufc_.form.num_coefficients(); ++i)
  {
    coefficients_data_size += ufc_.coefficient_elements[i]->space_dimension();
  }
  uint dofs_data_size = 0;
  for (uint i = 0; i < ufc_.form.rank(); ++i)
  {
    dofs_data_size += ufc_.local_dimensions[i];
  }

  // Separate real and uint data types to avoid copy of dof indices and casts
  // vertex coordinates + coefficients values
  r_packet_size_ = coordinates_data_size + coefficients_data_size;
  // local facet index + arguments dof indices
  u_packet_size_ = 1 + dofs_data_size;

  // Allocate data structures
  _set<uint> const& adj = distdata.get_adj_ranks();
  uint const num_shared_facets = distdata.num_shared();

  // Pack by adjacent rank
  uint offset = 0;
  for (_set<uint>::const_iterator it = adj.begin(); it != adj.end(); ++it)
  {
    rank_offsets_.insert(FacetOffsets(*it, offset));
    offset += distdata.shared_mapping().to(*it).size();
  }
  //
  dolfin_assert(offset == num_shared_facets);

  // Cache association of local facet index to offset in halo data
  _map<uint, uint> facet_offsets;
  for (SharedIterator sh(distdata); sh.valid(); ++sh)
  {
    uint const ark = *(distdata.get_shared_adj(sh.index()).begin());
    uint rank_offset = rank_offsets_[ark];

    // Maps local shared ordering to adjacent shared ordering
    Array<uint> const& adjmap = distdata.shared_mapping().to(ark);

    uint off0 = rank_offset + facet_offsets[ark];
    uint off1 = rank_offset + adjmap[facet_offsets[ark]];
    facet_map_[sh.index()] = FacetOffsets(off0, off1);

    // Increment current facet index for the given adjacent rank
    ++facet_offsets[ark];
  }
  //
  dolfin_assert(facet_map_.size() == num_shared_facets);

  r_data0_ = new real[num_shared_facets * r_packet_size_];
  u_data0_ = new uint[num_shared_facets * u_packet_size_];
  r_data1_ = new real[num_shared_facets * r_packet_size_];
  u_data1_ = new uint[num_shared_facets * u_packet_size_];

  // Fill data structures
  this->update(coefficients_, dof_map_set_);
}

//-----------------------------------------------------------------------------
void UFCHalo::update(Facet& facet)
{

  FacetMap::const_iterator it = facet_map_.find(facet.index());
  if (it == facet_map_.end())
  {
    if (facet.is_shared())
    {
      error("Shared facet index not found in halo data structure.");
    }
    else
    {
      error("Trying to fetch halo data for a non-shared facet.");
    }
  }

  // Update real data stored in halo data arrays
  real * r0 = &r_data0_[r_packet_size_ * it->second.first];
  real * r1 = &r_data1_[r_packet_size_ * it->second.second];

  // Update pointers to coordinates
  uint const gdim = mesh_.geometry_dimension();
  for (uint i = 0; i < mesh_.type().num_entities(0); ++i)
  {
    cell0.coordinates[i] = r0;
    r0 += gdim;
    cell1.coordinates[i] = r1;
    r1 += gdim;
  }

  // Update UFC expansion coefficients, needs copy for the moment
  for (uint i = 0; i < ufc_.form.num_coefficients(); ++i)
  {
    uint const spacedim = ufc_.coefficient_elements[i]->space_dimension();
    std::copy(r0, r0 + spacedim, macro_w[i]);
    r0 += spacedim;
    std::copy(r1, r1 + spacedim, macro_w[i] + spacedim);
    r1 += spacedim;
  }

  // Update uint data stored in halo data arrays:
  uint * u0 = &u_data0_[u_packet_size_ * it->second.first];
  uint * u1 = &u_data1_[u_packet_size_ * it->second.second];

  // Update local facet indices
  facet0 = *u0;
  ++u0;
  facet1 = *u1;
  ++u1;

  // Update UFC dofs indices for each dimension, needs copy for the moment
  for (uint i = 0; i < ufc_.form.rank(); ++i)
  {
    uint const localdim = ufc_.local_dimensions[i];
    std::copy(u0, u0 + localdim, macro_dofs[i]);
    u0 += localdim;
    std::copy(u1, u1 + localdim, macro_dofs[i] + localdim);
    u1 += localdim;
  }
}

//-----------------------------------------------------------------------------
void UFCHalo::update(Array<Coefficient*> const& coefficients,
                     DofMapSet const& dof_map_set)
{
  Mesh& mesh = dof_map_set[0].mesh();

  if (!mesh.is_distributed())
  {
    return;
  }

#ifdef DOLFIN_HAVE_MPI

  uint const tdim = mesh.topology_dimension();
  uint const gdim = mesh.geometry_dimension();
  uint const facet_dim = mesh.type().facet_dim();
  DistributedData& distdata = mesh.distdata()[facet_dim];

  // Exchange of data for contribution of halo macro elements

  if (coefficients.size() != ufc_.form.num_coefficients())
  {
    error("UFCHalo: invalid number of coefficients passed as arguments:\n"
          "Expected: %d ; Provided: %d",
          ufc_.form.num_coefficients(), coefficients.size());
  }

  //
  uint rank = dolfin::MPI::rank();
  uint pe_size = dolfin::MPI::size();

  // Loop over shared facets to collect data following shared iterator ordering
  uint const num_cell_vertices = mesh_.type().num_entities(0);
  for (FacetMap::const_iterator it = facet_map_.begin(); it != facet_map_.end();
      ++it)
  {
    Facet facet(mesh, it->first);

    // Create cell
    Cell cell(mesh, facet.entities(tdim)[0]);
    uint cell_facet = cell.index(facet);
    ufc_.cell.update(cell);

    // Set arrays offset
    real * r_entry = &r_data0_[it->second.first * r_packet_size_];
    uint * u_entry = &u_data0_[it->second.first * u_packet_size_];

    // Collect data for shared facet contribution
    // TODO: implement proper serialization functions
    for (uint i = 0; i < num_cell_vertices; ++i)
    {
      std::copy(ufc_.cell.coordinates[i], ufc_.cell.coordinates[i] + gdim,
                r_entry);
      r_entry += tdim;
    }

    // Interpolate coefficients on cell
    for (uint i = 0; i < ufc_.form.num_coefficients(); ++i)
    {
      coefficients[i]->interpolate(r_entry, ufc_.cell,
                                   *ufc_.coefficient_elements[i], cell,
                                   cell_facet);
      r_entry += ufc_.coefficient_elements[i]->space_dimension();
    }

    // Add local facet index
    *u_entry = cell_facet;
    ++u_entry;

    // Tabulate dofs for each dimension on macro element
    for (uint i = 0; i < ufc_.form.rank(); ++i)
    {
      dof_map_set[i].tabulate_dofs(u_entry, ufc_.cell, cell);
      u_entry += ufc_.local_dimensions[i];
    }
  }

  // Exchange data to fill halo data arrays: facet blocks are written directly
  for (int j = 1; j < (int) pe_size; ++j)
  {
    int src = (rank - j + pe_size) % pe_size;
    int dst = (rank + j) % pe_size;

    uint num_send_facets = distdata.shared_mapping().to(dst).size();
    uint num_recv_facets = distdata.shared_mapping().from(src).size();
    uint send_offset = rank_offsets_[dst];
    uint recv_offset = rank_offsets_[src];

    real * r_sendbuf = &r_data0_[send_offset * r_packet_size_];
    real * r_recvbuf = &r_data1_[recv_offset * r_packet_size_];
    int r_sendcount = num_send_facets * r_packet_size_;
    int r_recvcount = num_recv_facets * r_packet_size_;
    MPI::sendrecv( &r_sendbuf[0], r_sendcount, dst,
                   &r_recvbuf[0], r_recvcount, src, 1 );

    uint * u_sendbuf = &u_data0_[send_offset * u_packet_size_];
    uint * u_recvbuf = &u_data1_[recv_offset * u_packet_size_];
    int u_sendcount = num_send_facets * u_packet_size_;
    int u_recvcount = num_recv_facets * u_packet_size_;
    MPI::sendrecv( &u_sendbuf[0], u_sendcount, dst,
                   &u_recvbuf[0], u_recvcount, src, 1 );
  }

#else
  MAYBE_UNUSED(coefficients);
#endif
}

//-----------------------------------------------------------------------------
void UFCHalo::disp() const
{
  section("UFCHalo");
  prm("Facet map size", facet_map_.size());
  prm("Adjacent ranks", rank_offsets_.size());
  for (_map<uint, uint>::const_iterator it = rank_offsets_.begin();
       it != rank_offsets_.end(); ++it)
  {
    cout << "\tproc " << it->first << " : " << it->second << "\n";
  }
  prm("Size of real data packet", r_packet_size_);
  prm("Size of uint data packet", u_packet_size_);
  end();
}

//-----------------------------------------------------------------------------
void UFCHalo::clear()
{
  rank_offsets_.clear();
  facet_map_.clear();
  delete[] r_data0_;
  r_data0_ = nullptr;
  delete[] u_data0_;
  u_data0_ = nullptr;
  delete[] r_data1_;
  r_data1_ = nullptr;
  delete[] u_data1_;
  u_data1_ = nullptr;
}

} /* namespace dolfin */
