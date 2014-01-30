// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-27
// Last changed: 2014-01-27

#ifndef __UFL_TYPE_H_
#define __UFL_TYPE_H_

#include <dolfin/ufl/UFLObject.h>

#include <iostream>
#include <sstream>
#include <string>

namespace ufl
{

template<typename T>
class type : public Object
{

public:

  ///
  type<T>()
  {
  }

  /// Constructor with default representation for given type
  type<T>(T const& s) :
    val_(s),
    repr_(make_repr(s)),
    str_(make_str(s))
  {
  }

  ///
  ~type<T>()
  {
  }

  /// __repr__
  repr_t const repr() const;

  /// __str__
  std::string const str() const;

  ///
  virtual void display() const;

  ///
  repr_t const make_repr(
        std::vector<Object const *> const& prototype) const;

protected:

  ///
  repr_t const make_repr( T const& val) const;

  ///
  std::string const make_str( T const& val) const;

private:

  bool const is_string_type(std::string const& val) const;

  T val_;
  repr_t const repr_;
  std::string const str_;

};

//-----------------------------------------------------------------------------
template<typename T>
  Object::repr_t const type<T>::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
template<typename T>
  std::string const type<T>::str() const
{
  std::stringstream ss;
  ss << val_;
  return str_;
}

//-----------------------------------------------------------------------------
template<typename T>
  void type<T>::display() const
{
  std::cout << "Type" << std::endl;
  std::cout << ".str : " << (std::string) this->str() << std::endl;
  std::cout << ".repr: " << (std::string) this->repr() << std::endl;
  std::cout << std::endl;
}

//-----------------------------------------------------------------------------
template<typename T>
  Object::repr_t const type<T>::make_repr(
        std::vector<Object const *> const& prototype) const
{
  return Object::make_repr(prototype);
}

//-----------------------------------------------------------------------------
template<typename T>
  Object::repr_t const type<T>::make_repr(T const& val) const
{
  std::stringstream ret;
  ret << val;
  return ret.str();
}

//-----------------------------------------------------------------------------
template<typename T>
  std::string const type<T>::make_str( T const& val) const
{
  std::stringstream ret;
  ret << val;
  // Hack to handle string.
  std::string strval = ret.str();
  if(is_string_type(strval))
  {
    strval.erase(strval.begin());
    strval.erase(strval.end()-1);
  }
  return strval;
}

//-----------------------------------------------------------------------------
template<typename T>
  bool const type<T>::is_string_type(std::string const& val) const
{
  return (*val.begin() == '\'' && *val.rbegin() == '\'');
}


} /* namespace icorne */
#endif /* __UFL_TYPE_H_ */
