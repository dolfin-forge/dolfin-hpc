// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_CLASS_H
#define __DOLFIN_UFL_CLASS_H

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

  /// Return the name of the class, prefixed if needed
  std::string const name() const;

  /// Return the name of the class, prefixed
  std::pair<std::string, std::string> const prefixed_name() const;

  /// __repr__
  repr_t const& repr() const override = 0;

  /// __str__
  std::string const& str() const override = 0;

  /// __eq__
  bool operator ==(Class const& other) const;

  /// Create representation from a list of argument objects
  template<class OBJ>
    repr_t make_repr(std::vector<OBJ const *> const& args) const;

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
  explicit Class(std::string const& pre, std::string const& pos,
                 repr_t const& repr);

  ///
  ~Class() override;

  ///
  void display() const override;

  ///
  static std::string make_name(repr_t repr);

  ///
  repr_t const make_repr(Object const *& arg1) const;

  ///
  repr_t const make_repr(Object const *& arg1, Object const *& arg2) const;

  ///
  std::vector<repr_t> make_args_repr(repr_t const& repr, bool without_pre_pos =
                                         false) const override;

  ///
  repr_t const& arg(size_t i);

  ///
  std::vector<repr_t> const& args();

  ///
  void remove_pre_pos(repr_t const& repr, std::string& str, std::string& pre,
                      std::string& pos) const;

private:

  using CppProto = std::pair<std::pair<std::string, std::string>, std::vector<const Object *> >;

  ///
  CppProto make_proto(repr_t repr) const;

  //--- ATTRIBUTES ------------------------------------------------------------
  std::pair<std::string, std::string> const prefixed_name_;
  CppProto cpp_proto_;
  std::vector<repr_t> args_repr_;
  static repr_t const default_repr_;
  static std::string const default_str_;

};

//-----------------------------------------------------------------------------
inline bool Class::operator ==(Class const& other) const
{
  return (other.repr() == this->repr());
}

} /* namespace ufl */
#endif /* __DOLFIN_UFL_CLASS_H */
