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
  repr(Class const& owner, Object const& arg1);

  ///
  repr(Class const& owner, Object const& arg1, Object const& arg2);

  ///
  repr(Class const& owner, std::vector<Object const *> const& prototype);

  ///
  ~repr();

private:

  std::string const make_representation(
          Class const& owner, Object const& arg1);

  std::string const make_representation(
        Class const& owner, Object const& arg1, Object const& arg2);

  std::string const make_representation(
      Class const& owner, std::vector<Object const *> const& prototype);

};

} /* namespace icorne */
#endif /* __UFL_REPR_H_ */
