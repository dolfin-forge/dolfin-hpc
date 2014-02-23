// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_OBJECT_H_
#define __UFL_OBJECT_H_

#include <dolfin/ufl/UFLrepr.h>
#include <dolfin/common/types.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace ufl
{

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
  virtual repr_t const repr() const = 0;

  /// __str__
  virtual std::string const str() const = 0;

  ///
  virtual repr_t const make_repr(
      std::vector<Object const *> const& prototype) const = 0;

  /// __eq__
  virtual bool operator ==(Object const& other) const;

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
  virtual std::vector<repr_t> const make_args_repr(
      repr_t const& repr, bool const& without_pre_pos = false) const = 0;

  ///
  std::vector<Object const *> make_args(std::vector<repr_t> const& repr) const;

  /// Create from representation
  static Object * create(repr_t const& representation);

};

//-----------------------------------------------------------------------------
inline bool Object::operator ==(Object const& other) const
{
  return (other.repr() == this->repr());
}

//-----------------------------------------------------------------------------
inline void Object::display() const
{
  std::cout << "Object" << std::endl;
  std::cout << ".str : " << this->str() << std::endl;
  std::cout << ".repr: " << this->repr() << std::endl;
}

//-----------------------------------------------------------------------------
inline Object::repr_t const Object::make_repr(
    std::vector<Object const *> const& args) const
{
  std::cout << "Object::make_repr" << std::endl;
  std::stringstream ret;
  std::vector<Object const *>::const_iterator arg = args.begin();
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
inline std::vector<Object::repr_t> const Object::make_args_repr(
    repr_t const& repr, bool const& without_pre_pos) const
{
  //assumes repr to be a comma separated list
  std::vector<Object::repr_t> args;
  std::string str = repr;
  std::string delimiter = ", ";

  std::vector<char> open_delimiters;
  std::vector<char> close_delimiters;
  std::vector<size_t> open_delimiter_positions(3,0);
  std::vector<size_t> close_delimiter_positions(3,0);
  std::vector<size_t> n_open_delimiters(3,0);
  std::vector<size_t> n_close_delimiters(3,0);

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

  std::cout << "Parsing" << std::endl;
  while (scpos != MAX_STRING_LENGTH)
  {
    scpos = str.find(delimiter, currpos);
    std::string::iterator it = (scpos == std::string::npos ? str.end() : str.begin() + scpos) ;
    for(dolfin::uint i = 0; i<open_delimiters.size(); ++i)
    {
      open_delimiter_positions[i] = str.find(open_delimiters[i], currpos);
      close_delimiter_positions[i] = str.find(close_delimiters[i], currpos);
      n_open_delimiters[i] = std::count(str.begin(), it, open_delimiters[i]);
      n_close_delimiters[i] = std::count(str.begin(), it, close_delimiters[i]);
    }

    bool equal_number_open_close = true;
    for(dolfin::uint i = 0; i<n_open_delimiters.size(); ++i)
      if(n_open_delimiters[i] != n_close_delimiters[i])
        equal_number_open_close = false;

    open_pos = *std::min_element(open_delimiter_positions.begin(), open_delimiter_positions.end());
    dolfin::uint index = distance(open_delimiter_positions.begin(), std::find(open_delimiter_positions.begin(),
          open_delimiter_positions.end(), open_pos));
    close_pos = close_delimiter_positions[index];

    // Take into account of
    // - one argument
    // - class arguments
    // - tuple argument
    // - dict argument
    if (scpos == std::string::npos //no comma
        || (open_pos > scpos && equal_number_open_close) //comma is enclosed in a tuple or class
        || ((open_pos < scpos && close_pos < scpos) && equal_number_open_close)) //comma is after a tuple or a class
    {
      token = str.substr(0, scpos);

      if(scpos != std::string::npos)
        str.erase(0, scpos + delimiter.length());
      else
        str.erase(0, scpos);

      args.push_back(Object::repr_t(token));
      for(dolfin::uint i = 0; i<open_delimiters.size(); ++i)
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
//inline Object * Object::create(repr_t const& repr)
//{
//}

} /* namespace ufl */
#endif /* __UFL_OBJECT_H_ */
