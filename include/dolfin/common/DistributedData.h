//
//
//

#ifndef __DOLFIN_DISTRIBUTED_DATA_H
#define __DOLFIN_DISTRIBUTED_DATA_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>

namespace dolfin
{

/*
 *  @class  DistributedData
 */

struct DistributedData
{

  //--- Numbering
  bool valid_numbering;

  uint global_size;
  _map<uint, uint> global_indices;
  _map<uint, uint> local_indices;

  // Caching
  bool finalized;
  uint * cached_global_indices;

  //--- Ownership
  bool valid_ownership;
  Array<uint> ownership;

  // Shared
  _set<uint> shared;
  _set<uint> adjacent_ranks;
  _map<uint, _set<uint> > shared_adj;

  // Ghosts
  _set<uint> ghosts;
  _map<uint, uint> ghost_owner;

  //--- Adjacent mappings and reverse mappings
  bool valid_mapping;
  typedef _map<uint, std::pair< Array<uint>, Array<uint> > > AdjacentMapping;
  AdjacentMapping shared_mapping;
  AdjacentMapping ghost_mapping;

  //---------------------------------------------------------------------------
  inline void map(uint local_index, uint global_index)
  {
    global_indices.insert(std::pair<uint, uint>(local_index, global_index));
    local_indices.insert(std::pair<uint, uint>(global_index, local_index));
  }
  //---------------------------------------------------------------------------
  inline uint get_global(uint local_index)
  {
    return global_indices.find(local_index)->second;
  }
  //---------------------------------------------------------------------------
  inline uint get_local(uint global_index)
  {
    return global_indices.find(global_index)->second;
  }
  //---------------------------------------------------------------------------

};

} /* namespace dolfin */

#endif /* __DOLFIN_DISTRIBUTED_DATA_H */
