//
//
//

#ifndef __DOLFIN_LINEAR_DISTRIBUTION_H
#define __DOLFIN_LINEAR_DISTRIBUTION_H

#include <dolfin/common/types.h>
#include <dolfin/main/MPI.h>

#include <cmath>

namespace dolfin
{

struct LinearDistribution
{

  uint const global_size;
  uint const card;
  uint const rank;
  uint const L;
  uint const R;
  uint const offset;
  uint const size;

  ///
  LinearDistribution(uint global_size, uint card) :
      global_size(global_size),
      card(card),
      rank(MPI::processNumber()),
      L(std::floor((real) global_size / (real) card)),
      R(global_size % card),
      offset(rank * L + std::min(rank,R)),
      size((global_size + card - rank - 1) / card)
  {
  }

  ///
  inline bool in_range(uint index)
  {
    return ((offset <= index) && (index < offset + size));
  }

  ///
  inline uint owner(uint index) const
  {
    return static_cast<uint>(std::max(
        std::floor((real) index / (real) (L + 1)),
        std::floor((real) ((real) index - (real) R) / (real) L)));
  }

  ///
  void disp() const
  {
    section("LinearDistribution");
    message("global size    : %u", global_size);
    message("num partitions : %u", card);
    message("rank           : %u", rank);
    message("offset         : %u", offset);
    message("size           : %u", size);
    message("quotient size  : %u", L);
    message("remain         : %u", R);
    end();
  }

};

} /* namespace dolfin */

#endif /* __DOLFIN_LINEAR_DISTRIBUTION_H */
