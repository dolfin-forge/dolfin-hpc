// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_OBJECT_H
#define __DOLFIN_UFL_OBJECT_H

#include <dolfin/ufl/UFLrepr.h>

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

using dolfin::begin;
using dolfin::end;

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ufl
{

using dolfin::error;

/**
 *  DOCUMENTATION:
 *
 *  @class  Object
 *
 *  @brief  Provides an interface for Python objects from UFL.
 */

class Object
{
  static std::size_t const MAX_STRING_LENGTH = 2048;

public:

  typedef ufl::repr repr_t;

  /// __repr__
  virtual repr_t const& repr() const = 0;

  /// __str__
  virtual std::string const& str() const = 0;

  ///
  template<class OBJ>
    repr_t make_repr(std::vector<OBJ const *> const& prototype) const;

  /// __eq__
  bool operator ==(Object const& other) const;

  /// __neq__
  bool operator !=(Object const& other) const;

protected:

  ///
  Object()
  {
  }

  ///
  virtual ~Object()
  {
  }

  ///
  virtual void display() const = 0;

  ///
  virtual std::vector<repr_t> make_args_repr(
      repr_t const& repr, bool without_pre_pos = false) const = 0;

  ///
  std::vector<Object const *> make_args(std::vector<repr_t> const& repr) const;

  /// Create from representation
  static Object const * create(repr_t const& representation);
};

//-----------------------------------------------------------------------------
inline bool Object::operator ==(Object const& other) const
{
  return (other.repr() == this->repr());
}

//-----------------------------------------------------------------------------
inline bool Object::operator !=(Object const& other) const
{
  return (other.repr() != this->repr());
}

//-----------------------------------------------------------------------------
inline void Object::display() const
{
  begin("UFL::Object");
  std::cout << "Object" << std::endl;
  std::cout << ".str : " << this->str() << std::endl;
  std::cout << ".repr: " << this->repr() << std::endl;
}

//-----------------------------------------------------------------------------
template<class OBJ>
  inline Object::repr_t Object::make_repr(
      std::vector<OBJ const *> const& args) const
  {
    std::stringstream ret;
    if (args.size() == 0) return ret.str();
    typename std::vector<OBJ const *>::const_iterator arg = args.begin();
    ret << (*arg)->repr();
    for (++arg; arg != args.end(); ++arg)
    {
      ret << ", " << (*arg)->repr();
    }
    return ret.str();
  }

//-----------------------------------------------------------------------------
inline std::vector<Object const *> Object::make_args(
    std::vector<repr_t> const& repr) const
{
  std::vector<Object const *> args;
  for (std::vector<repr_t>::const_iterator it = repr.begin(); it != repr.end();
      ++it)
  {
    args.push_back(this);
  }
  return args;
}

//-----------------------------------------------------------------------------
inline std::vector<Object::repr_t> Object::make_args_repr(
    repr_t const& repr, bool) const
{
  //assumes repr to be a comma separated list
  std::vector<Object::repr_t> args;
  std::string str = repr;

  size_t star_pos = 0;
  while (star_pos != std::string::npos && star_pos < str.length())
  {
    star_pos = str.find("*");
    if (star_pos != std::string::npos) str.erase(star_pos, 1);
  }

  std::string delimiter = ", ";

  std::vector<char> open_delimiters;
  std::vector<char> close_delimiters;
  std::vector<size_t> open_delimiter_positions(3, 0);
  std::vector<size_t> close_delimiter_positions(3, 0);
  std::vector<size_t> n_open_delimiters(3, 0);
  std::vector<size_t> n_close_delimiters(3, 0);

  open_delimiters.push_back('(');
  open_delimiters.push_back('[');
  open_delimiters.push_back('{');

  close_delimiters.push_back(')');
  close_delimiters.push_back(']');
  close_delimiters.push_back('}');

  size_t scpos = 0;
  size_t currpos = 0;
  size_t open_pos = 0;
  size_t close_pos = 0;

  std::string token;
  while (scpos != std::string::npos || scpos < str.size())
  {
    scpos = str.find(delimiter, currpos);
    std::string::iterator it = (
        scpos == std::string::npos ? str.end() : str.begin() + scpos);
    for (dolfin::uint i = 0; i < open_delimiters.size(); ++i)
    {
      open_delimiter_positions[i] = str.find(open_delimiters[i], currpos);
      close_delimiter_positions[i] = str.find(close_delimiters[i], currpos);
#ifdef __SUNPRO_CC
      n_open_delimiters[i] = 0;
      n_close_delimiters[i] = 0;
      std::count(str.begin(), it, open_delimiters[i], n_open_delimiters[i]);
      std::count(str.begin(), it, close_delimiters[i], n_close_delimiters[i]);
#else
      n_open_delimiters[i] = std::count(str.begin(), it, open_delimiters[i]);
      n_close_delimiters[i] = std::count(str.begin(), it, close_delimiters[i]);
#endif
    }

    bool equal_number_open_close = true;
    for (dolfin::uint i = 0; i < n_open_delimiters.size(); ++i)
      if (n_open_delimiters[i] != n_close_delimiters[i]) equal_number_open_close =
          false;

    open_pos = *std::min_element(open_delimiter_positions.begin(),
                                 open_delimiter_positions.end());
#ifdef __SUNPRO_CC
    dolfin::uint index = 0;
    std::distance(open_delimiter_positions.begin(),
        std::find(open_delimiter_positions.begin(),
            open_delimiter_positions.end(), open_pos), index);
#else
    dolfin::uint index = std::distance(
        open_delimiter_positions.begin(),
        std::find(open_delimiter_positions.begin(),
                  open_delimiter_positions.end(), open_pos));
#endif
    close_pos = close_delimiter_positions[index];

    // Take into account of
    // - one argument
    // - class arguments
    // - tuple argument
    // - dict argument
    if (scpos == std::string::npos  //no comma
    || (open_pos > scpos && equal_number_open_close)  //comma is enclosed in a tuple or class
        || ((open_pos < scpos && close_pos < scpos) && equal_number_open_close))  //comma is after a tuple or a class
    {
      token = str.substr(0, scpos);

      if (scpos != std::string::npos) str.erase(0, scpos + delimiter.length());
      else str.erase(0, scpos);

      args.push_back(Object::repr_t(token));
      for (dolfin::uint i = 0; i < open_delimiters.size(); ++i)
      {
        open_delimiter_positions[i] = str.find(open_delimiters[i], currpos);
        close_delimiter_positions[i] = str.find(close_delimiters[i], currpos);
      }
      currpos = 0;
    }
    else
    {
      currpos = scpos + 1;
    }
  }
  return args;
}

//-----------------------------------------------------------------------------
inline Object const* Object::create(repr_t const&)
{
  //FIXME complete this list
//  std::string name = Class::make_name(repr);
//  if (name == "Sum")
//  {
//    return new Sum(repr);
//  }
//  else if (name == "Product")
//  {
//    return new Product(repr);
//  }
//  else if (name == "Division")
//  {
//    return new Division(repr);
//  }
//  else if (name == "Power")
//  {
//    return new Power(repr);
//  }
//  else if (name == "Abs")
//  {
//    return new Abs(repr);
//  }
//  else if (name == "Argument")
//  {
//    return new Argument(repr);
//  }
//  else if (name == "Coefficient")
//  {
//    return new Coefficient(repr);
//  }
//  else if (name == "Constant")
//  {
//    return new Constant(repr);
//  }
//  else if (name == "VectorConstant")
//  {
//    return new VectorConstant(repr);
//  }
//  else if (name == "TensorConstant")
//  {
//    return new TensorConstant(repr);
//  }
//  else if (name == "CoefficientDerivative")
//  {
//    return new CoefficientDerivative(repr);
//  }
//  else if (name == "SpatialDerivative")
//  {
//    return new SpatialDerivative(repr);
//  }
//  else if (name == "VariableDerivative")
//  {
//    return new VariableDerivative(repr);
//  }
//  else if (name == "Grad")
//  {
//    return new Grad(repr);
//  }
//  else if (name == "Div")
//  {
//    return new Div(repr);
//  }
//  else if (name == "NablaGrad")
//  {
//    return new NablaGrad(repr);
//  }
//  else if (name == "NablaDiv")
//  {
//    return new NablaDiv(repr);
//  }
//  else if (name == "Curl")
//  {
//    return new Curl(repr);
//  }
//  else if (name == "Index")
//  {
//    return new Index(repr);
//  }
//  else if (name == "FixedIndex")
//  {
//    return new FixedIndex(repr);
//  }
//  else if (name == "MultiIndex")
//  {
//    return new MultiIndex(repr);
//  }
//  else if (name == "IndexSum")
//  {
//    return new IndexSum(repr);
//  }
//  else if (name == "Indexed")
//  {
//    return new Indexed(repr);
//  }
//  else if (name == "ComponentTensor")
//  {
//    return new ComponentTensor(repr);
//  }
//  else if (name == "Tuple")
//  {
//    return new Tuple(repr);
//  }
//  else if (name == "Label")
//  {
//    return new Label(repr);
//  }
//  else if (name == "Variable")
//  {
//    return new Variable(repr);
//  }
//  else
//  {
//    error("Unknown type of ufl::Object: '" + name + "'");
//  }

  return NULL;
}
} /* namespace ufl */
#endif /* __DOLFIN_UFL_OBJECT_H */
