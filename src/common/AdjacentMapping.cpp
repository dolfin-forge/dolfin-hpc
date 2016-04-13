//
//
//

#include <dolfin/common/AdjacentMapping.h>

#include <dolfin/main/MPI.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
SharedMapping::SharedMapping(DistributedData const& data) :
    data_(data),
    mappings_(),
    send_min_(0),
    send_max_(0)
{
  if(!data_.is_finalized())
  {
    error("SharedMapping : distributed data is not finalized");
  }

#if HAVE_MPI

  uint const rank = dolfin::MPI::processNumber();
  uint const pe_size = dolfin::MPI::numProcesses();

  //
  for (SharedIterator it(data); !it.end(); ++it)
  {
    _set<uint> const& adjs = it.adj();
    for (_set<uint>::const_iterator a = adjs.begin(); a != adjs.end(); ++a)
    {
      mappings_[*a].send.push_back(it.global_index());
    }
  }
  dolfin_assert(mappings_.size() == data.get_adj_ranks().size());
  send_max_ = 0;
  send_min_ = data.num_shared();
  for (_map<uint, AdjacentMapping>::iterator it = mappings_.begin();
       it != mappings_.end(); ++it)
  {
    dolfin_assert(it->first != rank);
    // Update bounds
    send_max_ = std::max(send_max_, (uint) it->second.send.size());
    send_min_ = std::min(send_min_, (uint) it->second.send.size());
    // Resize buffer
    it->second.recv.resize(it->second.send.size());
  }

  // Exchange lists of global indices
  MPI_Status status;
  uint dst;
  uint src;
  int recv_count;
  for (int j = 1; j < (int) pe_size; ++j)
  {
    src = (rank - j + pe_size) % pe_size;
    dst = (rank + j) % pe_size;

    uint * senddst = NULL;
    uint sendcnt = 0;
    _map<uint, AdjacentMapping>::iterator adjdst = mappings_.find(dst);
    if (adjdst != mappings_.end())
    {
      senddst = &adjdst->second.send[0];
      sendcnt = adjdst->second.send.size();
    }

    uint * recvbuf = NULL;
    uint recvcnt = 0;
    _map<uint, AdjacentMapping>::iterator adjsrc = mappings_.find(src);
    if (adjsrc != mappings_.end())
    {
      recvbuf = &adjsrc->second.recv[0];
      recvcnt = adjsrc->second.recv.size();
    }

    MPI_Sendrecv(&senddst[0], sendcnt, MPI_UNSIGNED, dst, 1, &recvbuf[0],
                 recvcnt, MPI_UNSIGNED, src, 1, dolfin::MPI::DOLFIN_COMM,
                 &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);
  }

  // Mappings contain global indices, map back to local
  message("Map back to local");
  for (_map<uint, AdjacentMapping>::iterator it = mappings_.begin();
       it != mappings_.end(); ++it)
  {
    data_.get_local(it->second.send.size(), &it->second.send[0], &it->second.send[0]);
    data_.get_local(it->second.recv.size(), &it->second.recv[0], &it->second.recv[0]);
  }

#endif /* HAVE_MPI */

}
//-----------------------------------------------------------------------------
SharedMapping::SharedMapping(SharedMapping const& other) :
    data_(other.data_),
    mappings_(other.mappings_),
    send_min_(other.send_min_),
    send_max_(other.send_max_)
{
}
//-----------------------------------------------------------------------------
SharedMapping::~SharedMapping()
{
}
//-----------------------------------------------------------------------------
SharedMapping& SharedMapping::operator=(SharedMapping const& other)
{
  // Do not allow assignment
  return *this;
}
//-----------------------------------------------------------------------------
Array<uint> const& SharedMapping::to(uint rank) const
{
  _map<uint, AdjacentMapping>::const_iterator it = mappings_.find(rank);
  if (it == mappings_.end())
  {
    error("SharedMapping : invalid adjacent %u", rank);
  }
  return it->second.send;
}
//-----------------------------------------------------------------------------
Array<uint> const& SharedMapping::from(uint rank) const
{
  _map<uint, AdjacentMapping>::const_iterator it = mappings_.find(rank);
  if (it == mappings_.end())
  {
    error("SharedMapping : invalid adjacent %u", rank);
  }
  return it->second.recv;
}
//-----------------------------------------------------------------------------
void SharedMapping::disp() const
{
  section("SharedMapping");
  message("number of adjacents : %u", mappings_.size());
  message("minimum size        : %u", send_min_);
  message("maximum size        : %u", send_max_);
  end();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
