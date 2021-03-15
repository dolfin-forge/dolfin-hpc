// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_SUB_FUNCTION_H
#define __DOLFIN_SUB_FUNCTION_H

#include <dolfin/common/types.h>
#include <dolfin/function/Function.h>
#include <dolfin/log/log.h>

namespace dolfin
{

class Function;

//-----------------------------------------------------------------------------

/**
 *  @class  SubFunction
 *
 *  @brief  This class represents a sub function of a (discrete function) to
 *          enable expressions like:
 *
 *            Function w;
 *            Function u = w[0];
 *            Function p = w[1];
 *
 *          without needing to create and destroy temporaries.
 *          No data is created until a Function is assigned to a SubFunction, at
 *          which point the data needed to represent the subfunction is created.
 */

class SubFunction
{

public:
  /// Create sub function
  SubFunction( Function & f, size_t i )
    : f_( &f )
    , i_( i )
  {
  }

  /// Destructor
  ~SubFunction() = default;

  /// Return global function
  auto function() const -> Function &;

  /// Return index of the sub function
  auto index() const -> size_t;

  /// Display basic information
  void disp() const;

private:
  /// Create empty sub function
  SubFunction()

  {
  }

  // Pointer to discrete function
  Function * const f_ { nullptr };

  // Sub function index
  size_t const i_ { 0 };
};

//-----------------------------------------------------------------------------

inline auto SubFunction::function() const -> Function &
{
  return *f_;
}

//-----------------------------------------------------------------------------

inline auto SubFunction::index() const -> size_t
{
  return i_;
}

//-----------------------------------------------------------------------------

inline void SubFunction::disp() const
{
  section( "SubFunction" );
  prm( "Index", this->index() );
  function().disp();
  end();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_SUB_FUNCTION_H */
