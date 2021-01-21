
#ifndef __DOLFIN_ADJACENT_MAPPING_H
#define __DOLFIN_ADJACENT_MAPPING_H

#include <dolfin/common/Array.h>
#include <dolfin/common/DistributedData.h>

namespace dolfin
{

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
  SharedMapping( DistributedData const & data );

  ///
  SharedMapping( SharedMapping const & other ) = default;

  ///
  ~SharedMapping() = default;

  /// Do not allow assignment
  auto operator=( SharedMapping const & other ) -> SharedMapping & = delete;

  ///
  auto to( uint rank ) const -> Array< uint > const &;

  ///
  auto from( uint rank ) const -> Array< uint > const &;

  ///
  void disp() const;

private:
  struct AdjacentMapping
  {
    Array< uint > send;
    Array< uint > recv;

    ///
    AdjacentMapping()
      : send()
      , recv()
    {
    }
  };

private:

  DistributedData const &       data_;
  _map< uint, AdjacentMapping > mappings_;
  uint                          send_min_;
  uint                          send_max_;
};

} /* namespace dolfin */

#endif /* __DOLFIN_ADJACENT_MAPPING_H */
