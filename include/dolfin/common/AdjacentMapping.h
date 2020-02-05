
#ifndef __DOLFIN_ADJACENT_MAPPING_H
#define __DOLFIN_ADJACENT_MAPPING_H

#include <dolfin/common/Array.h>
#include <dolfin/common/DistributedData.h>

namespace dolfin
{

/**
 *  @class  AdjacentMapping
 *
 *  @brief
 *
 */

struct AdjacentMapping
{
  Array<uint> send;
  Array<uint> recv;

  ///
  AdjacentMapping() :
      send(),
      recv()
  {
  }

};

/**
 *  @class  SharedMapping
 *
 *  @brief
 *
 */

class SharedMapping
{

public:

  ///
  SharedMapping(DistributedData const& data);

  ///
  SharedMapping(SharedMapping const& other);

  ///
  ~SharedMapping();

  ///
  Array<uint> const& to(uint rank) const;

  ///
  Array<uint> const& from(uint rank) const;

  ///
  void disp() const;

private:

  /// Do not allow assignment
  SharedMapping& operator=(SharedMapping const& other);

  DistributedData const& data_;
  _map<uint, AdjacentMapping> mappings_;
  uint send_min_;
  uint send_max_;

};

} /* namespace dolfin */

#endif /* __DOLFIN_ADJACENT_MAPPING_H */
