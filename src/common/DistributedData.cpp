// Copyright (C) 2016 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//

#include <dolfin/common/DistributedData.h>

#include <dolfin/log/log.h>
#include <dolfin/main/MPI.h>

#include <cstring>

namespace dolfin
{

//-----------------------------------------------------------------------------
DistributedData::DistributedData() :
    valid_numbering(false),
    valid_ownership(false),
    valid_adjacency(false),
    rank_(MPI::processNumber()),
    pe_size_(MPI::numProcesses()),
    offset_(0),
    range_size_(0),
    global_size_(0),
    finalized_(false),
    global_(),
    local_(),
    adjacents_(),
    shared_(),
    ghost_(),
    cache_size_(0),
    cached_numbering_(NULL),
    cached_ownership_(NULL)
{
}
//-----------------------------------------------------------------------------
DistributedData::DistributedData(DistributedData const& other)
{
  *this = other;
}
//-----------------------------------------------------------------------------
DistributedData::~DistributedData()
{
  clear();
}
//-----------------------------------------------------------------------------
DistributedData& DistributedData::operator=(DistributedData const& other)
{
  if (this != &other)
  {
    clear();

    valid_numbering = other.valid_numbering;
    valid_ownership = other.valid_ownership;

    rank_ = other.rank_;
    pe_size_ = other.pe_size_;
    offset_ = other.offset_;
    range_size_ = other.range_size_;
    global_size_ = other.global_size_;

    finalized_ = other.finalized_;

    if (finalized_)
    {
      // Avoid copy of global-to-local unordered map
      cache_size_ = other.cache_size_;
      dolfin_assert(other.cache_size_ == other.local_.size());
      if (cache_size_ >  0)
      {
        dolfin_assert(other.cached_numbering_ !=  NULL);
        dolfin_assert(other.cached_ownership_ !=  NULL);
        cached_numbering_ = new uint[cache_size_];
        cached_ownership_ = new uint[cache_size_];
        for (uint i = 0; i < other.local_.size(); ++i)
        {
          cached_numbering_[i] = other.cached_numbering_[i];
          local_[cached_numbering_[i]] = i;
          cached_ownership_[i] = other.cached_ownership_[i];
        }
        dolfin_assert(local_.size() == other.local_.size());
      }
    }
    else
    {
      cache_size_ = other.cache_size_;
      if (other.cached_numbering_ != NULL)
      {
        dolfin_assert(cache_size_ > 0);
        cached_numbering_ = new uint[cache_size_];
        std::copy(other.cached_numbering_,
                  other.cached_numbering_ + cache_size_, cached_numbering_);
      }
      if (other.cached_ownership_ != NULL)
      {
        dolfin_assert(cache_size_ > 0);
        cached_ownership_ = new uint[cache_size_];
        std::copy(other.cached_ownership_,
                  other.cached_ownership_ + cache_size_, cached_ownership_);
      }
      global_ = other.global_;
      local_ = other.local_;
    }

    adjacents_ = other.adjacents_;
    shared_ = other.shared_;
    ghost_ = other.ghost_;
  }
  return *this;
}
//-----------------------------------------------------------------------------
bool DistributedData::operator==(DistributedData const& other) const
{
  return true;
}
//-----------------------------------------------------------------------------
bool DistributedData::operator!=(DistributedData const& other) const
{
  return !(*this == other);
}
//-----------------------------------------------------------------------------
void DistributedData::clear()
{
  delete [] cached_ownership_;
  delete [] cached_numbering_;
  cache_size_ = 0;
  ghost_.clear();
  shared_.clear();
  adjacents_.clear();
  local_.clear();
  global_.clear();
  finalized_ = false;
  global_size_ = 0;
  range_size_ = 0;
  offset_ = 0;
  pe_size_ = MPI::numProcesses();
  rank_ = MPI::processNumber();
  //
  valid_numbering = false;
  valid_ownership = false;
  valid_adjacency = false;
}
//-----------------------------------------------------------------------------
void DistributedData::finalize()
{
  if (finalized_)
  {
    if (local_.size() > 0)
    {
      if (cached_numbering_ == NULL)
      {
        error("DistributedData : data is finalized but empty numbering cache");
      }
      if (cached_ownership_ == NULL)
      {
        error("DistributedData : data is finalized but empty ownership cache");
      }
    }
    // Be conservative for the moment
    error("DistributedData : trying to finalize data twice");
  }
  else
  {
    message(1, "DistributedData : finalize");

    // Check consistency
    if (shared_.size() > local_.size())
    {
      error("DistributedData : shared size is greater than mapping size");
    }
    if (ghost_.size() > shared_.size())
    {
      error("DistributedData : ghost size is greater than shared size");
    }

    // Cache provided mapping for numbering and ownership
    if (global_.size() > 0)
    {
      // Generate caching for existing mapping
      if (cached_numbering_ != NULL)
      {
        error("DistributedData : data not finalized but numbering is cached");
      }
      if (cached_ownership_ != NULL)
      {
        error("DistributedData : data not finalized but ownership is cached");
      }

      message(1, "DistributedData : cache mapping of size %u", global_.size());
      if (global_.size() != local_.size())
      {
        error("DistributedData : size mismatch between index mappings %u != %u",
              global_.size(), local_.size());
      }

      // Cache numbering and ownership
      cache_size_ = global_.size();
      cached_numbering_ = new uint[cache_size_];
      cached_ownership_ = new uint[cache_size_];
      for (_map<uint, uint>::iterator it = global_.begin();
           it != global_.end(); ++it)
      {
        cached_numbering_[it->first] = it->second;
        cached_ownership_[it->first] = pe_size_;
      }
      global_.clear();

      // Update ownership for shared entities
      for (_map<uint, _set<uint> >::const_iterator it = shared_.begin();
           it != shared_.end(); ++it)
      {
        cached_ownership_[it->first] = rank_;
      }

      // Update ownership for ghost entities
      for (_map<uint, uint>::const_iterator it = ghost_.begin();
           it != ghost_.end(); ++it)
      {
        cached_ownership_[it->first] = it->second;
      }
    }

    // No mapping provided but a process range can be set.
    // Generate a linear mapping from process range: entities are not ghosted.
    if ((local_.size() == 0) && (range_size_ > 0 || cache_size_ > 0))
    {
      if (shared_.size() > 0)
      {
        error("DistributedData : no mapping but entities set as shared");
      }

      // Either local size or the range can be used
      if (cache_size_ > 0)
      {
        if ((range_size_ > 0) && (range_size_ != cache_size_))
        {
          error("DistributedData : size mismatch between local size and range");
        }
        range_size_ = cache_size_;
      }
      else
      {
        cache_size_ = range_size_;
        cached_numbering_ = new uint[range_size_];
        cached_ownership_ = new uint[range_size_];
      }

      // Numbering incrementally and set all as owned
      MPI::processOffset(range_size_, offset_);
      uint global = offset_;
      for (uint local = 0; local < range_size_; ++local, ++global)
      {
        cached_numbering_[local] = global;
        local_[global] = local;
      }
      std::fill_n(cached_ownership_, range_size_, pe_size_);
      // Global renumbering is not necessary
      valid_numbering = true;

      message(1, "DistributedData : generated linear mapping in range [%u, %u[",
              offset_, offset_+range_size_);
    }

    // At this point mappings exist and are cached.
    // For the sake of completeness let us check the consistency
    if (local_.size() != cache_size_)
    {
      error("DistributedData : size mismatch between local-to-global (%u) and "
            "and global-to-local (%u) mappings", cache_size_, local_.size());
    }

    // Set data range if needed
    uint owned_size = local_.size() - ghost_.size();
    if (range_size_ == 0)
    {
      range_size_ = owned_size;
      MPI::processOffset(range_size_, offset_);
    }
    else if (range_size_ != owned_size)
    {
      error("DistributedData : data range does not match provided range");
    }

    // Check global size anyway
    uint range_sum;
    MPI::numGlobalSum(range_size_, range_sum);
    // Check that computed value matches the former value such that the sum of
    // ranges is indeed equal to the previously set global size
    if (global_size_ == 0)
    {
      global_size_ = range_sum;
    }
    else if (global_size_ != range_sum)
    {
      error("DistributedData : sum of range not equal to global size");
    }

    //
    finalized_ = true;
  }
}
//-----------------------------------------------------------------------------
bool DistributedData::empty() const
{
  return (local_.size() == 0 && global_.size() == 0 && shared_.size() == 0
          && ghost_.size() == 0);
}
//-----------------------------------------------------------------------------
uint DistributedData::capacity() const
{
  return local_.size();
}
//-----------------------------------------------------------------------------
bool DistributedData::is_finalized() const
{
  return finalized_;
}
//-----------------------------------------------------------------------------
uint DistributedData::offset() const
{
  return offset_;
}
//-----------------------------------------------------------------------------
uint DistributedData::range_size() const
{
  return range_size_;
}
//-----------------------------------------------------------------------------
bool DistributedData::in_range(uint global_index) const
{
  dolfin_assert(global_size_ > 0);
  dolfin_assert(offset_+ range_size_ <= global_size_);
  return (offset_ <= global_index && global_index < offset_+ range_size_);
}
//-----------------------------------------------------------------------------
uint DistributedData::local_size() const
{
  return local_.size();
}
//-----------------------------------------------------------------------------
uint DistributedData::global_size() const
{
  if (global_size_ < local_.size())
  {
    error("DistributedData : global size has not been set or is invalid");
  }
  return global_size_;
}
//-----------------------------------------------------------------------------
void DistributedData::set_range(uint num_owned, uint num_global /* = 0 */ )
{
  if (finalized_)
  {
    error("DistributedData : setting range to a finalized data");
  }
  // Check if there exists a mapping already
  if (!this->empty())
  {
    error("DistributedData : setting range to a non-empty data");
  }
  // Check if any existing global size matches the provided value
  if ((range_size_ > 0) && (range_size_ != num_owned))
  {
    error("DistributedData : setting different range than existing");
  }
  // Check if the provided range size is consistent with local size
  if (num_owned < cache_size_)
  {
    error("DistributedData : provided range is smaller than local size ");
  }
  // Set range in any case
  range_size_ = num_owned;
  MPI::processOffset(range_size_, offset_);
  // Set the global size if provided otherwise compute it
  if (num_global > 0)
  {
    // Check if any existing global size matches the provided value
    if ((global_size_ > 0) && (global_size_ != num_global))
    {
      error("DistributedData : setting different global size than existing");
    }
    global_size_ = num_global;
  }
  else
  {
    uint range_sum;
    MPI::numGlobalSum(range_size_, range_sum);
    // Check that computed value matches the former value such that the sum of
    // ranges is indeed equal to the previously set global size
    if ((global_size_ > 0) && (global_size_ != range_sum))
    {
      error("DistributedData : sum of range not equal to global size");
    }
  }
}
//-----------------------------------------------------------------------------
void DistributedData::set_size(uint num_local, uint num_global /* = 0 */ )
{
  if (finalized_)
  {
    error("DistributedData : setting size to a finalized data");
  }
  // Do not allow setting the size of of a non-empty distributed data
  if (!this->empty())
  {
    error("DistributedData : setting size to a non-empty data");
  }
  // Do not allow resetting the size of data
  if ((cached_numbering_ != NULL) && (cache_size_ != num_local))
  {
    error("DistributedData : setting different size to data");
  }
  // If range is known check consistency, range is zero by default
  if (num_local < range_size_)
  {
    error("DistributedData : setting smaller local size than existing range");
  }
  // Set the global size if provided
  if (num_global > 0)
  {
    if ((global_size_ > 0) && (global_size_ != num_global))
    {
      error("DistributedData : setting different global size than existing");
    }
    global_size_ = num_global;
  }
  // Use caching already as the data size is known and save on map lookups
  cache_size_ = num_local;
  // Create arrays if they do not exist
  if ((cached_numbering_ == NULL) & (cache_size_ > 0))
  {
    cached_numbering_ = new uint[cache_size_];
    // Set to an undefined value
    std::fill_n(cached_numbering_, cache_size_, DOLFIN_UINT_UNDEF);
  }
  if ((cached_ownership_ == NULL) & (num_local > 0))
  {
    cached_ownership_ = new uint[cache_size_];
    // Set initial ownership to owner
    std::fill_n(cached_ownership_, cache_size_, pe_size_);
  }
}
//-----------------------------------------------------------------------------
uint DistributedData::has_local(uint local_index) const
{
  if (cached_numbering_ != NULL)
  {
    dolfin_assert(local_index < cache_size_);
    return cached_numbering_[local_index] != DOLFIN_UINT_UNDEF;
  }
  return (global_.count(local_index) > 0);
}
//-----------------------------------------------------------------------------
uint DistributedData::get_global(uint local_index) const
{
  if (cached_numbering_ != NULL)
  {
    dolfin_assert(local_index < cache_size_);
    dolfin_assert(cached_numbering_[local_index] != DOLFIN_UINT_UNDEF);
    return cached_numbering_[local_index];
  }
  dolfin_assert(global_.count(local_index) > 0);
  return global_.find(local_index)->second;
}
//-----------------------------------------------------------------------------
uint DistributedData::has_global(uint global_index) const
{
  return (local_.count(global_index) > 0);
}
//-----------------------------------------------------------------------------
uint DistributedData::get_local(uint global_index) const
{
  dolfin_assert(local_.count(global_index) > 0);
  return local_.find(global_index)->second;
}
//---------------------------------------------------------------------------
void DistributedData::set_map(uint local_index, uint global_index)
{
  dolfin_assert(!finalized_);
  dolfin_assert(local_.count(global_index) == 0);
  local_.insert(std::pair<uint, uint>(global_index, local_index));
  if (cached_numbering_ != NULL)
  {
    dolfin_assert(local_index < cache_size_);
    cached_numbering_[local_index] = global_index;
  }
  else
  {
    dolfin_assert(global_.count(local_index) == 0);
    global_.insert(std::pair<uint, uint>(local_index, global_index));
  }
}
//-----------------------------------------------------------------------------
void DistributedData::remap_numbering(Array<uint> const& mapping)
{
  if (!finalized_)
  {
    error("DistributedData : re-mapping numbering requires finalized data");
  }

  if (mapping.size() != MPI::numProcesses())
  {
    error("DistributedData : numbering re-mapping array has invalid size");
  }

  // Update numbering
  dolfin_assert(global_.size() == 0);
  dolfin_assert(cached_numbering_ != NULL);
  for (_map<uint, uint>::iterator it = local_.begin(); it != local_.end();
       ++it)
  {
    cached_numbering_[it->second] = mapping[it->second];
    it->second = mapping[it->second];
  }

  // Update shared entities
  dolfin_assert(cached_ownership_ != NULL);
  _map<uint, _set<uint> > shared;
  for (_map<uint, _set<uint> >::const_iterator it = shared_.begin();
       it != shared_.end(); ++it)
  {
    shared[mapping[it->first]] = it->second;
    cached_ownership_[mapping[it->first]] = rank_;
  }
  shared_ = shared;

  // Update ghost entities
  dolfin_assert(cached_ownership_ != NULL);
  _map<uint, uint> ghost;
  for (_map<uint, uint>::const_iterator it = ghost_.begin(); it != ghost_.end();
       ++it)
  {
    ghost[mapping[it->first]] = it->second;
    cached_ownership_[mapping[it->first]] = it->second;
  }
  ghost_ = ghost;
}
//-----------------------------------------------------------------------------
void DistributedData::renumber_global()
{
  if (!finalized_)
  {
    error("DistributedData : global renumbering requires finalized data");
  }

#if HAVE_MPI

  message(1, "DistributedData : renumber global");

  _map<uint, uint> local_mapping;
  Array<uint> * sendbuf = new Array<uint> [pe_size_];

  // Re-index owned entities and collect ghosted entities per owner
  dolfin_assert(!(local_.size() > 0 && cached_numbering_ == NULL));
  dolfin_assert(!(local_.size() > 0 && cached_ownership_ == NULL));
  uint index = offset_;
  for (uint i = 0; i < local_.size(); ++i)
  {
    dolfin_assert(cached_ownership_[i] <= pe_size_);
    if (cached_ownership_[i] == pe_size_)
    {
      cached_numbering_[i] = index;
      local_mapping[index] = i;
      ++index;
    }
    else
    {
      sendbuf[cached_ownership_[i]].push_back(cached_numbering_[i]);
    }
  }

  // Exchange data and set numbering
  MPI_Status status;
  uint src;
  uint dst;

  // Maximum number of received entities is the number of owned shared but the
  // issue is that some previous were written in a way that does not ensure
  // appropriate setting of owned shared entities
  uint recvsize = sendbuf[0].size();
  for (uint i = 1; i < pe_size_; ++i)
  {
    recvsize = std::max(recvsize, (uint) sendbuf[i].size());
  }
  MPI::numGlobalMax(recvsize, recvsize);
  uint * recvbuf = (recvsize == 0 ? NULL : new uint[recvsize]);
  uint * sendbck = (recvsize == 0 ? NULL : new uint[recvsize]);
  uint const num_ghost = ghost_.size();
  uint * recvbck = (num_ghost == 0 ? NULL : new uint[num_ghost]);

  int recvcount;
  for (uint j = 1; j < pe_size_; ++j)
  {
    src = (rank_ - j + pe_size_) % pe_size_;
    dst = (rank_ + j) % pe_size_;

    MPI_Sendrecv(&sendbuf[dst][0], sendbuf[dst].size(), MPI_UNSIGNED, dst, 1,
                 &recvbuf[0], recvsize, MPI_UNSIGNED, src, 1,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recvcount);

    for (int k = 0; k < recvcount; ++k)
    {
      dolfin_assert(local_.count(recvbuf[k]) > 0);
      uint const local_index = local_.find(recvbuf[k])->second;
      sendbck[k] = cached_numbering_[local_index];

      // Entity is not marked as shared: invalidate ownership
      if (cached_ownership_[local_index] == pe_size_)
      {
        valid_ownership = false;
      }
      else if (cached_ownership_[local_index] != rank_)
      {
        error("DistributedData : entity %u should be owned but is set ghost");
      }
    }

    MPI_Sendrecv(&sendbck[0], recvcount, MPI_UNSIGNED, src, 2,
                 &recvbck[0], sendbuf[dst].size(), MPI_UNSIGNED, dst, 2,
                 MPI::DOLFIN_COMM, &status);

    for (uint k = 0; k < sendbuf[dst].size(); ++k)
    {
      uint const local_index = local_.find(sendbuf[dst][k])->second;
      dolfin_assert(local_mapping.count(recvbck[k]) == 0);
      cached_numbering_[local_index] = recvbck[k];
      local_mapping[recvbck[k]] = local_index;
    }
  }

  delete[] recvbck;
  delete[] sendbck;
  delete[] recvbuf;
  delete[] sendbuf;

  local_ = local_mapping;

  message("Done");
#endif /* HAVE_MPI */

}
//-----------------------------------------------------------------------------
bool DistributedData::has_adj_rank(uint rank) const
{
  return (adjacents_.count(rank) > 0);
}
//-----------------------------------------------------------------------------
uint DistributedData::num_adj_ranks() const
{
  return adjacents_.size();
}
//-----------------------------------------------------------------------------
_set<uint> const& DistributedData::get_adj_ranks() const
{
  return adjacents_;
}
//-----------------------------------------------------------------------------
uint DistributedData::get_owner(uint local_index) const
{
  dolfin_assert(this->has_local(local_index));
  if (cached_ownership_ != NULL)
  {
    return (cached_ownership_[local_index] == pe_size_ ?
              rank_ : cached_ownership_[local_index]);
  }
  _map<uint, uint>::const_iterator it = ghost_.find(local_index);
  if (it == ghost_.end())
  {
    return rank_;
  }
  return it->second;
}
//-----------------------------------------------------------------------------
bool DistributedData::is_owned(uint local_index) const
{
  dolfin_assert(this->has_local(local_index));
  if (cached_ownership_ != NULL)
  {
    return (cached_ownership_[local_index] == pe_size_ ||
            cached_ownership_[local_index] == rank_);
  }
  return (ghost_.count(local_index) == 0);
}
//-----------------------------------------------------------------------------
bool DistributedData::is_shared(uint local_index) const
{
  dolfin_assert(this->has_local(local_index));
  if (cached_ownership_ != NULL)
  {
    return (cached_ownership_[local_index] < pe_size_);
  }
  return (shared_.count(local_index) == 0);
}
//-----------------------------------------------------------------------------
bool DistributedData::is_ghost(uint local_index) const
{
  dolfin_assert(this->has_local(local_index));
  if (cached_ownership_ != NULL)
  {
    return (cached_ownership_[local_index] < pe_size_ &&
            cached_ownership_[local_index] != rank_);
  }
  return (ghost_.count(local_index) > 0);
}
//-----------------------------------------------------------------------------
uint DistributedData::num_owned() const
{
  return (local_.size() - ghost_.size());
}
//-----------------------------------------------------------------------------
uint DistributedData::num_shared() const
{
  return (shared_.size());
}
//-----------------------------------------------------------------------------
uint DistributedData::num_ghost() const
{
  return (ghost_.size());
}
//-----------------------------------------------------------------------------
void DistributedData::remap_ownership(Array<uint> const& mapping)
{
  dolfin_assert(finalized_);
  if (mapping.size() != MPI::numProcesses())
  {
    error("DistributedData : ownership re-mapping array has invalid size");
  }

  // Update current rank
  rank_ = mapping[rank_];

  // Update adjacent ranks
  _set<uint> adjs;
  for (_set<uint>::const_iterator it = adjacents_.begin();
       it != adjacents_.end(); ++it)
  {
    adjs.insert(mapping[*it]);
  }
  adjacents_ = adjs;

  // Update shared entities ownership
  for (_map<uint, _set<uint> >::iterator it = shared_.begin();
       it != shared_.end(); ++it)
  {
    // Update adjacent
    _set<uint> adj;
    for (_set<uint>::const_iterator a = it->second.begin();
         a != it->second.end(); ++a)
    {
      adj.insert(mapping[*a]);
    }
    it->second = adj;

    // Update cached owner
    if (cached_ownership_ != NULL)
    {
      cached_ownership_[it->first] = rank_;
    }
  }

  // Update ghost entities ownership
  for (_map<uint, uint>::iterator it = ghost_.begin(); it != ghost_.end();
       ++it)
  {
    // Update owner
    it->second = mapping[it->second];

    // Update cached owner
    if (cached_ownership_ != NULL)
    {
      cached_ownership_[it->first] = mapping[it->second];
    }
  }
}
//-----------------------------------------------------------------------------
_set<uint> const& DistributedData::get_shared_adj(uint local_index) const
{
  dolfin_assert(shared_.count(local_index) > 0);
  return shared_.find(local_index)->second;
}
//-----------------------------------------------------------------------------
void DistributedData::get_common_adj(uint n, uint const indices[],
                                     _set<uint>& adjs) const
{
  if (n == 0)
  {
    adjs.clear();
    return;
  }
  dolfin_assert(shared_.count(indices[0]) > 0);
  adjs = shared_.find(indices[0])->second;
  for (uint i = 1; (adjs.size()> 0)&&(i < n); ++i)
  {
    dolfin_assert(shared_.count(indices[i]) > 0);
    _set<uint> const& adjx = shared_.find(indices[i])->second;
    for(_set<uint>::iterator it = adjs.begin(); it != adjs.end();)
    {
      if (adjx.count(*it) == 0)
      {
        adjs.erase(it++);
      }
      else
      {
        ++it;
      }
    }
  }
}
//-----------------------------------------------------------------------------
void DistributedData::set_shared(uint local_index)
{
  dolfin_assert(!finalized_);
  dolfin_assert(this->has_local(local_index));
  dolfin_assert(shared_.count(local_index) == 0);
  if ((cached_ownership_ != NULL)
      && (cached_ownership_[local_index] == pe_size_))
  {
    cached_ownership_[local_index] = rank_;
  }
  // As explained we allow setting an entity as shared without adjacent only if
  // the entity is not shared already.
  dolfin_assert(global_.count(local_index) > 0);
  if (shared_[local_index].size() > 0)
  {
    error("DistributedData : cannot set_shared on entities with adjacents");
  }
}
//-----------------------------------------------------------------------------
void DistributedData::set_shared_adj(uint local_index, uint adj)
{
  dolfin_assert(!finalized_);
  dolfin_assert(this->has_local(local_index));
  dolfin_assert(adj != rank_);
  dolfin_assert(adj < pe_size_);
  if ((cached_ownership_ != NULL)
      && (cached_ownership_[local_index] == pe_size_))
  {
    cached_ownership_[local_index] = rank_;
  }
  shared_[local_index].insert(adj);
  adjacents_.insert(adj);
}
//-----------------------------------------------------------------------------
void DistributedData::setall_shared_adj(uint local_index, _set<uint> const& adjs)
{
  dolfin_assert(!finalized_);
  dolfin_assert(this->has_local(local_index));
  dolfin_assert(adjs.count(rank_) == 0);
  if ((cached_ownership_ != NULL)
      && (cached_ownership_[local_index] == pe_size_))
  {
    cached_ownership_[local_index] = rank_;
  }
  shared_[local_index] = adjs;
  adjacents_.insert(adjs.begin(), adjs.end());
}
//-----------------------------------------------------------------------------
void DistributedData::set_ghost(uint local_index, uint owner)
{
  dolfin_assert(!finalized_);
  dolfin_assert(this->has_local(local_index));
  dolfin_assert(owner != rank_);
  dolfin_assert(owner < pe_size_);
  if (cached_ownership_ != NULL)
  {
    cached_ownership_[local_index] = owner;
  }
  shared_[local_index].insert(owner);
  ghost_[local_index] = owner;
  adjacents_.insert(owner);
}
//-----------------------------------------------------------------------------
void DistributedData::disp() const
{
  section("DistributedData");
  message("rank        : %8u", rank_);
  message("offset      : %8u", offset_);
  message("range size  : %8u", range_size_);
  message("local  size : %8u", local_.size());
  message("global size : %8u", global_size_);
  message("finalized   : %8u", finalized_);
  /*
  if (finalized_)
  {
    section("entities    :");
    for (uint i = 0; i < this->local_size(); ++i)
    {
      message("%8u -> %8u -> %8u : %8u", i, cached_numbering_[i],
              local_.find(cached_numbering_[i])->second, this->get_owner(i));
    }
    end();
  }
  */
  end();
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
