// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-27
// Last changed: 2014-01-27

#ifndef __UFL_REPR_H_
#define __UFL_REPR_H_

#include <string>
#include <vector>

namespace ufl
{

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

  ///
  ~repr();

  ///
  operator const char *() { return this->c_str(); }

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
#endif /* __UFL_REPR_H_ */
