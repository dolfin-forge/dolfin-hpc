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
    repr_(type::make_repr(s)),
    str_(type::make_str(s))
  {
  }

  /// Constructor with default representation for given type
  type<T>(Object::repr_t const& repr) :
    val_(type::make_val(repr)),
    repr_(repr),
    str_(type::make_str(val_))
  {
  }

  ///
  ~type<T>()
  {
  }

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  virtual void display() const;

  ///
  repr_t make_repr(std::vector<Object const *> const& prototype) const;

  ///
//  repr_t const make_repr(
//        std::vector<Expression const *> const& prototype) const;

  operator T() const { return val_;}

protected:

  ///
  repr_t make_repr( T const& val) const;

  ///
  T make_val(repr_t const& repr) const;

  ///
  std::string make_str( T const& val) const;

  ///
  std::vector<Object::repr_t> make_args_repr(repr_t const& repr,
                                             bool without_pre_pos = false) const;

private:

  bool is_string_type(std::string const& val) const;

  T val_;
  repr_t const repr_;
  std::string const str_;

};

//-----------------------------------------------------------------------------
template<typename T>
  Object::repr_t const& type<T>::repr() const
{
  return repr_;
}

//-----------------------------------------------------------------------------
template<typename T>
  std::string const& type<T>::str() const
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
  Object::repr_t type<T>::make_repr(
        std::vector<Object const *> const& prototype) const
{
  return Object::make_repr(prototype);
}

//-----------------------------------------------------------------------------
template<typename T>
  Object::repr_t type<T>::make_repr(T const& val) const
{
  std::stringstream ret;
  ret << val;
  return ret.str();
}

//-----------------------------------------------------------------------------
template<typename T>
  T type<T>::make_val(repr_t const& repr) const
{
  T val;
  std::stringstream ss;
  ss << repr;
  ss >> std::ws >> val; // Important to not ignore whitespaces
  return val;
}

//-----------------------------------------------------------------------------
template<typename T>
  std::string type<T>::make_str( T const& val) const
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
  std::vector<Object::repr_t> type<T>::make_args_repr(repr_t const& repr,
                                                      bool without_pre_pos) const
{
  return std::vector<Object::repr_t>();
}

//-----------------------------------------------------------------------------
template<typename T> bool type<T>::is_string_type(std::string const& val) const
{
  return (*val.begin() == '\'' && *val.rbegin() == '\'');
}


} /* namespace ufl */
#endif /* __UFL_TYPE_H_ */
