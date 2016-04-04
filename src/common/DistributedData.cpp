//
//
//

#include <dolfin/common/DistributedData.h>

#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>
#include <dolfin/main/MPI.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------
DistributedData::DistributedData() :
    valid_numbering(false),
    valid_ownership(false),
    rank_(MPI::processNumber()),
    pe_size_(MPI::numProcesses()),
    offset_(0),
    global_size_(0),
    finalized_(false),
    global_indices_(),
    local_indices_(),
    cached_global_indices_(NULL),
    ownership_(),
    adjacents_(),
    shared_(),
    shared_adj_(),
    ghosts_(),
    ghost_owner_()
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
  delete [] cached_global_indices_;
}
//-----------------------------------------------------------------------------
DistributedData& DistributedData::operator=(DistributedData const& other)
{
  if (this != &other)
  {
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
}
//-----------------------------------------------------------------------------
void DistributedData::finalize()
{
  if (finalized_)
  {
    if (cached_global_indices_ == NULL)
    {
      error("DistributedData : data is finalized but cache is empty");
    }
  }
  else
  {
    if (cached_global_indices_ != NULL)
    {
      error("DistributedData : data is not finalized and cache is not empty");
    }
    if (global_indices_.size() != local_indices_.size())
    {
      error("DistributedData : size mismatch between index mappings %u != %u",
            global_indices_.size(), local_indices_.size());
    }
    cached_global_indices_ = new uint[global_indices_.size()];
    for (_map<uint, uint>::iterator it = global_indices_.begin();
         it != global_indices_.end(); ++it)
    {
      cached_global_indices_[it->first] = it->second;
    }
    global_indices_.clear();
    //
    finalized_ = true;
  }
}
//-----------------------------------------------------------------------------
uint DistributedData::capacity() const
{
  return ownership_.size();
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
  return (offset_ <= global_index && global_index < offset_+ range_size_);
}
//-----------------------------------------------------------------------------
uint DistributedData::local_size() const
{
  return local_indices_.size();
}
//-----------------------------------------------------------------------------
uint DistributedData::global_size() const
{
  return global_size_;
}
//---------------------------------------------------------------------------
void DistributedData::set_map(uint local_index, uint global_index)
{
  dolfin_assert(!finalized_);
  dolfin_assert(global_indices_.count(local_index) == 0);
  global_indices_.insert(std::pair<uint, uint>(local_index, global_index));
  dolfin_assert(local_indices_.count(global_index) == 0);
  local_indices_.insert(std::pair<uint, uint>(global_index, local_index));
}
//-----------------------------------------------------------------------------
void DistributedData::remap(Array<uint> const& mapping)
{
  dolfin_assert(finalized_);
  dolfin_assert(mapping.size() == MPI::numProcesses());

  // Update current rank
  rank_ = mapping[rank_];
  // Update shared entities ownership
  dolfin_assert(mapping.size() == MPI::numProcesses());
  for (_map<uint, _set<uint> >::iterator it = shared_adj_.begin();
       it != shared_adj_.end(); ++it)
  {
    // Update owner
    ownership_[it->first] = mapping[ownership_[it->first]];
    // Update adjacent
    _set<uint> adj;
    for(_set<uint>::const_iterator adjit = it->second.begin();
        adjit != it->second.end(); ++adjit)
    {
      adj.insert(mapping[*adjit]);
    }
    it->second = adj;
  }
}
//-----------------------------------------------------------------------------
void DistributedData::renumber_global()
{
  dolfin_assert(finalized_);
  error("DistributedData::renumber_global TBD");
}
//-----------------------------------------------------------------------------
uint DistributedData::has_local(uint local_index) const
{
  if (cached_global_indices_ != NULL)
  {
    return (local_index < local_indices_.size());
  }
  return (global_indices_.count(local_index) > 0);
}
//-----------------------------------------------------------------------------
uint DistributedData::get_global(uint local_index) const
{
  dolfin_assert(this->has_local(local_index));
  return (cached_global_indices_ != NULL ?
            cached_global_indices_[local_index] :
            global_indices_.find(local_index)->second);
}
//-----------------------------------------------------------------------------
uint DistributedData::has_global(uint global_index) const
{
  return (local_indices_.count(global_index) > 0);
}
//-----------------------------------------------------------------------------
uint DistributedData::get_local(uint global_index) const
{
  dolfin_assert(local_indices_.count(global_index) > 0);
  return local_indices_.find(global_index)->second;
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
  return (ownership_[local_index] == pe_size_ ? rank_
                                              : ownership_[local_index]);
}
//-----------------------------------------------------------------------------
bool DistributedData::is_owned(uint local_index) const
{
  dolfin_assert(ownership_.size() == global_indices_.size());
  dolfin_assert(local_index < ownership_.size());
  return (ownership_[local_index] == pe_size_ ||
          ownership_[local_index] == rank_);
}
//-----------------------------------------------------------------------------
bool DistributedData::is_shared(uint local_index) const
{
  dolfin_assert(ownership_.size() == global_indices_.size());
  dolfin_assert(local_index < ownership_.size());
  return (ownership_[local_index] < pe_size_);
}
//-----------------------------------------------------------------------------
bool DistributedData::is_ghost(uint local_index) const
{
  dolfin_assert(ownership_.size() == global_indices_.size());
  dolfin_assert(local_index < ownership_.size());
  return (ownership_[local_index] < pe_size_ &&
          ownership_[local_index] != rank_);
}
//-----------------------------------------------------------------------------
uint DistributedData::num_owned() const
{
  return (global_indices_.size() - ghosts_.size());
}
//-----------------------------------------------------------------------------
uint DistributedData::num_shared() const
{
  return (shared_.size());
}
//-----------------------------------------------------------------------------
uint DistributedData::num_ghost() const
{
  return (ghosts_.size());
}
//-----------------------------------------------------------------------------
_set<uint> const& DistributedData::get_shared_adj(uint local_index) const
{
  dolfin_assert(is_shared(local_index));
  dolfin_assert(shared_adj_.count(local_index) > 0);
  return shared_adj_.find(local_index)->second;
}
//-----------------------------------------------------------------------------
void DistributedData::set_shared(uint local_index)
{
  dolfin_assert(!finalized_);
  dolfin_assert(this->has_local(local_index));
  // Only set shared if the entity is not
  if(ownership_[local_index] == pe_size_)
  {
    shared_.insert(local_index);
    ownership_[local_index] = rank_;
  }
  dolfin_assert(this->is_shared(local_index));
}
//-----------------------------------------------------------------------------
void DistributedData::set_shared_adj(uint local_index, uint adj)
{
  dolfin_assert(!finalized_);
  dolfin_assert(this->has_local(local_index));
  dolfin_assert(this->is_shared(local_index));
  shared_adj_[local_index].insert(adj);
  adjacents_.insert(adj);
}
//-----------------------------------------------------------------------------
void DistributedData::setall_shared_adj(uint local_index, _set<uint> const& adjs)
{
  dolfin_assert(!finalized_);
  dolfin_assert(this->has_local(local_index));
  dolfin_assert(adjs.count(rank_) == 0);
  dolfin_assert(this->is_shared(local_index));
  shared_adj_[local_index].insert(adjs.begin(), adjs.end());
  adjacents_.insert(adjs.begin(), adjs.end());
}
//-----------------------------------------------------------------------------
void DistributedData::set_ghost(uint local_index, uint owner)
{
  dolfin_assert(!finalized_);
  dolfin_assert(this->has_local(local_index));
  dolfin_assert(owner != rank_);
  if(ownership_[local_index] == pe_size_)
  {
    shared_.insert(local_index);
    ghosts_.insert(local_index);
  }
  else if (ownership_[local_index] == rank_)
  {
    ghosts_.insert(local_index);
  }
  shared_adj_[local_index].insert(owner);
  adjacents_.insert(owner);
  ownership_[local_index] = owner;
  dolfin_assert(this->is_ghost(local_index));
}
//-----------------------------------------------------------------------------
void DistributedData::disp() const
{
  section("DistributedData");
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
