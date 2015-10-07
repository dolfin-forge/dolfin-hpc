// Copyright (C) 2013 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-05-24 (merged from branch larcher)
// Last changed: 2013-05-24

#ifndef __DOLFIN_EXPRESSION_H
#define __DOLFIN_EXPRESSION_H

#include <dolfin/common/types.h>
#include <dolfin/mesh/SubDomain.h>

namespace dolfin
{

/// Derived classes should implement the copy constructor if the allocate
/// new objects.

class Expression
{
public:

  /// Create user-defined function
  Expression()
  {
  }

  /// Return the rank of the value space
  virtual uint rank() const = 0;

  /// Return the dimension of the value space for axis i
  virtual uint dim(uint i) const = 0;

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const = 0;

protected:

  /// Destructor
  virtual ~Expression()
  {
  }

private:

};

class ScalarExpression : public Expression
{
public:

  /// Create user-defined function
  ScalarExpression() :
      Expression()
  {
  }

  /// Destructor
  ~ScalarExpression()
  {
  }

  /// Return the rank of the value space
  uint rank() const
  {
    return 0;
  }

  /// Return the dimension of the value space for axis i
  uint dim(uint i) const
  {
    return ( i > 0 ? 0 : 1 );
  }

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const = 0;

};

class VectorExpression : public Expression
{
public:

  /// Create user-defined function
  VectorExpression(uint value_dim) :
      Expression(),
      value_dim_(value_dim)
  {
  }

  /// Destructor
  ~VectorExpression()
  {
  }

  /// Return the rank of the value space
  uint rank() const
  {
    return 1;
  }

  /// Return the dimension of the value space for axis i
  uint dim(uint i) const
  {
    switch (i)
      {
      case 0:
        return value_dim_;
        break;
      case 1:
        return 1;
        break;
      default:
        return 0;
        break;
      }
    return value_dim_;
  }

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const = 0;

private:

  uint const value_dim_;

};

class TensorExpression : public Expression
{
public:

  /// Create user-defined function for 2nd order tensor
  TensorExpression(uint dim0, uint dim1) :
      Expression(),
      value_shape_(2)
  {
    value_shape_[0] = dim0;
    value_shape_[1] = dim1;
  }

  /// Create user-defined function for general tensor
  TensorExpression(Array<uint> value_shape) :
      Expression(),
      value_shape_(value_shape)
  {
  }

  /// Destructor
  ~TensorExpression()
  {
  }

  /// Return the rank of the value space
  uint rank() const
  {
    return value_shape_.size();
  }

  /// Return the dimension of the value space for axis i
  uint dim(uint i) const
  {
#ifdef __sgi
    dolfin_assert(i < value_shape_.size());
    return value_shape_[i];
#else
    return value_shape_.at(i);
#endif
  }

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const = 0;

private:

  Array<uint> value_shape_;

};

class RealReference : public ScalarExpression
{
public:

  /// Create user-defined function
  RealReference(real const& r) :
      ScalarExpression(),
      r_(r)
  {
  }

  /// Destructor
  ~RealReference()
  {
  }

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const
  {
    values[0] = r_;
  }

private:
  real const& r_;

};

class IndicatorExpression : public ScalarExpression
{
public:

  IndicatorExpression(SubDomain const& sub_domain, real value) :
      ScalarExpression(),
      sd_(sub_domain),
      value_(value)
  {
  }

  /// Destructor
  ~IndicatorExpression()
  {
  }

  inline void eval(real* values, const real* x) const
  {
    if (sd_.inside(x, true))
    {
      values[0] = value_;
    }
    else
    {
      values[0] = 0.0;
    }
  }

private:

  SubDomain const& sd_;
  real value_;

};

}

#endif
