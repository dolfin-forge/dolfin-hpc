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

  // Collect entities by adjacent rank
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

  //
  uint const rank = dolfin::MPI::processNumber();
  uint const pe_size = dolfin::MPI::numProcesses();
  MPI_Request * request = new MPI_Request[mappings_.size()];
  MPI_Status status;
  int recvcount;
  uint i = 0;
  for (_map<uint, AdjacentMapping>::iterator it = mappings_.begin();
       it != mappings_.end(); ++it, ++i)
  {
    dolfin_assert(it->first != rank);
    // Update bounds
    send_max_ = std::max(send_max_, (uint) it->second.send.size());
    send_min_ = std::min(send_min_, (uint) it->second.send.size());
    //
    MPI_Isend(&it->second.send[0], it->second.send.size(), MPI_UNSIGNED,
              it->first, 0, MPI::DOLFIN_COMM, &request[i]);
    // Resize buffer
    it->second.recv.resize(it->second.send.size());
    MPI_Irecv(&it->second.recv[0], it->second.recv.size(), MPI_UNSIGNED,
              it->first, 0, MPI::DOLFIN_COMM, &request[i]);
  }
  for (_map<uint, AdjacentMapping>::iterator it = mappings_.begin();
       it != mappings_.end(); ++it, ++i)
  {
    data_.get_local(it->second.send.size(), &it->second.send[0],
                    &it->second.send[0]);
  }
  i = 0;
  for (_map<uint, AdjacentMapping>::iterator it = mappings_.begin();
       it != mappings_.end(); ++it, ++i)
  {
    MPI_Wait(&request[i],&status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recvcount);
    if(recvcount != it->second.recv.size())
    {
      error("AdjacentMapping : inconsistent count %u from rank %u: expected %u",
            rank, recvcount, it->second.recv.size());
    }
    data_.get_local(it->second.recv.size(), &it->second.recv[0],
                    &it->second.recv[0]);
  }
  delete [] request;

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
