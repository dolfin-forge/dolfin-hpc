// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_CLASS_H_
#define __UFL_CLASS_H_

#include <string>
#include <vector>

#include <dolfin/ufl/UFLObject.h>

#include <dolfin/common/types.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Class
 *
 *  @brief  Provides an interface for Python classes from UFL.
 */

class Class : public Object
{

public:

  ///
  std::string const& name() const;

  /// __repr__
  virtual repr_t const repr() const = 0;

  /// __str__
  virtual std::string const str() const = 0;

  /// __eq__
  virtual bool operator ==(Class const& other) const;

  ///
  repr_t const make_repr(std::vector<Object const *> const& args) const;

protected:

  ///
  Class();

  ///
  Class(std::string const& name);
  
  ///
  Class(std::string const& pre, std::string const& pos);

  ///
  explicit Class(std::string const& name, repr_t const& repr);

  ///
  explicit Class(std::string const& pre, std::string const& pos, repr_t const& repr);

  ///
  virtual ~Class();

  ///
  virtual void display() const;

  ///
  static std::string const make_name(repr_t repr);

  ///
  repr_t const make_repr(Object const *& arg1) const;

  ///
  repr_t const make_repr(Object const *& arg1, Object const *& arg2) const;

  ///
  std::vector<repr_t> const make_args_repr(repr_t const& repr) const;

  ///
  repr_t const& arg(size_t i);

  ///
  std::vector<repr_t> const& args();

private:

  typedef std::pair<std::pair<std::string const, std::string const> const, std::vector<Object const*> > CppProto;

  ///
  CppProto make_proto(repr_t repr) const;

  //--- ATTRIBUTES ------------------------------------------------------------
  CppProto cpp_proto_;
  std::vector<repr_t> args_repr_;
  std::string const pre_;
  std::string const pos_;
  static repr_t const default_repr_;
  static std::string const default_str_;

};

//-----------------------------------------------------------------------------
inline bool Class::operator ==(Class const& other) const
{
  return (other.repr() == this->repr());
}

//-----------------------------------------------------------------------------
class ValueArray : public std::vector<dolfin::uint>
{

public:

  ///
  ValueArray();

  ///
  ValueArray(dolfin::uint const i);

  ///
  ValueArray(dolfin::uint const k, dolfin::uint const i);

  ///
  ~ValueArray();

  ///
  std::string const str() const;

private:

};

} /* namespace ufl */
#endif /* __UFL_CLASS_H */
