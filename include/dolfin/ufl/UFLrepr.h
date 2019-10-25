// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_REPR_H
#define __DOLFIN_UFL_REPR_H

#include <string>
#include <vector>

namespace ufl
{

//class Expression;
class Class;
class Object;

class repr : public std::string
{

public:

  ///
  repr();

  ///
  repr(std::string const& s);

  ///
  repr(const char *& s);

  ///
  repr(Class const& owner, Object const& arg1);

  ///
  repr(Class const& owner, Object const& arg1, Object const& arg2);

  ///
  repr(Class const& owner, Object const& arg1, Object const& arg2,
       Object const& arg3);

  ///
  repr(Class const& owner, Object const& arg1, Object const& arg2,
       Object const& arg3, Object const& arg4);

  ///
  repr(Class const& owner, Object const& arg1, Object const& arg2,
       Object const& arg3, Object const& arg4, Object const& arg5);

  ///
  repr(Class const& owner, std::vector<Object const *> const& prototype);

  //FIXME: Is it really necessary to add the same function also for Expressions?
//  repr(Expression const& owner, std::vector<Expression const *> const& prototype);

  template<class OBJ>
    repr(Class const& owner, std::vector<OBJ const *> const& prototype);

  ///
  ~repr();

  ///
  operator const char *()
  {
    return this->c_str();
  }

private:

  static std::string const make_repr(Class const& owner, Object const& arg1);

  static std::string const make_repr(Class const& owner, Object const& arg1,
                                     Object const& arg2);

  static std::string const make_repr(Class const& owner, Object const& arg1,
                                     Object const& arg2, Object const& arg3);

  static std::string const make_repr(Class const& owner, Object const& arg1,
                                     Object const& arg2, Object const& arg3,
                                     Object const& arg4);

  static std::string const make_repr(Class const& owner, Object const& arg1,
                                     Object const& arg2, Object const& arg3,
                                     Object const& arg4, Object const& arg5);
};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_REPR_H */
