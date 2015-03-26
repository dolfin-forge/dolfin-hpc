// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
//
// First added:  2015-03-12
// Last changed: 2015-03-12

#ifndef __CONSTANT_H
#define __CONSTANT_H

#include "Function.h"

namespace dolfin
{

class Constant : public Function
{

public:

  /// Create constant real function on given mesh
  Constant(Mesh& mesh, real val);

  ///
  ~Constant();
  
  //--- INTERFACE -------------------------------------------------------------

  /// Return the rank of the value space
  uint rank() const;

  /// Return the dimension of the value space for axis i
  uint dim(uint i) const;

  /// Evaluate function at given point
  void eval(real* values, const real* x) const;

  //---------------------------------------------------------------------------
  
  /// Assign constant real value (to all the components)
  Constant& operator= (real const& val);

private:

  Array<real> value_;

};

//-----------------------------------------------------------------------------
inline Constant::Constant(Mesh& mesh, real val) :
    Function(mesh),
    value_(static_cast<uint>(1), static_cast<real>(val))
{
}

//-----------------------------------------------------------------------------
inline Constant::~Constant()
{
}

//-----------------------------------------------------------------------------
inline void Constant::eval(real* values, const real* x) const
{
  values[0] = value_[0];
}

//-----------------------------------------------------------------------------
inline uint Constant::rank() const
{
  return 0;
}

//-----------------------------------------------------------------------------
inline uint Constant::dim(uint i) const
{
  return 1;
}

//-----------------------------------------------------------------------------
inline Constant& Constant::operator= (real const& val)
{
  std::fill(value_.begin(), value_.end(), val);
  return *this;
}

}

#endif /* __CONSTANT_H */
