// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_CLASS_H_
#define __UFL_CLASS_H_

#include <string>

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLCell
 *
 *  @brief  Provides an interface for Python objects from UFL.
 */

class UFLClass
{

public:

  /// __repr__
  virtual std::string const repr() const = 0;

  /// __str__
  virtual std::string const str() const = 0;

protected:

  UFLClass();

  virtual ~UFLClass();

private:

  static std::string const default_repr_;
  static std::string const default_str_;

};

} /* namespace dolfin */
#endif /* __UFL_CLASS_H */
