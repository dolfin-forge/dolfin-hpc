// Copyright (C) 2013 Aurélien Larcher
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-05-24 (merged from branch larcher)
// Last changed: 2013-05-24

#ifndef __EXPRESSION_H
#define __EXPRESSION_H

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

class ScalarExpression: public Expression
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
    return 1;
  }

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const = 0;

};

class VectorExpression: public Expression
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
    return value_dim_;
  }

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const = 0;

private:

  uint const value_dim_;

};

class RealReference: public Expression
{
public:

  /// Create user-defined function
  RealReference(real const& r) :
      Expression(),
      r_(r)
  {
  }

  /// Destructor
  ~RealReference()
  {
  }

  /// Return the rank of the value space
  virtual uint rank() const
  {
    return 0;
  }

  /// Return the dimension of the value space for axis i
  virtual uint dim(uint i) const
  {
    return 1;
  }

  /// Evaluate function at given point
  virtual void eval(real* values, const real* x) const
  {
    values[0] = r_;
  }

private:
  real const& r_;

};

class IndicatorExpression: public Expression
{
public:

  IndicatorExpression(SubDomain const& sub_domain, real value) :
      Expression(),
      sd_(sub_domain),
      value_(value)
  {
  }

  /// Destructor
  ~IndicatorExpression()
  {
  }

  inline uint rank() const
  {
    return 0;
  }

  inline uint dim(uint i) const
  {
    return 1;
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
