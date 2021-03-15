
#ifndef __DOLFIN_ADJACENT_MAPPING_H
#define __DOLFIN_ADJACENT_MAPPING_H

#include <dolfin/common/DistributedData.h>

#include <vector>

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
  auto to( size_t rank ) const -> std::vector< size_t > const &;

  ///
  auto from( size_t rank ) const -> std::vector< size_t > const &;

  ///
  void disp() const;

private:
  struct AdjacentMapping
  {
    std::vector< size_t > send;
    std::vector< size_t > recv;

    ///
    AdjacentMapping()
      : send()
      , recv()
    {
    }
  };

private:
  DistributedData const &         data_;
  _map< size_t, AdjacentMapping > mappings_;
  size_t                          send_min_;
  size_t                          send_max_;
};

} /* namespace dolfin */

#endif /* __DOLFIN_ADJACENT_MAPPING_H */
