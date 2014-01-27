// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-27
// Last changed: 2014-01-27

#ifndef __UFL_TYPE_H_
#define __UFL_TYPE_H_

#include <string>

namespace ufl
{

template<typename T>
class Type : public Object
{

public:

  ///
  Type<T>()
  {
  }

  ///
  Type<T>(T const& s) :
    val_(s)
  {
  }

  ///
  ~Type<T>()
  {
  }

  /// __repr__
  repr_t const repr() const;

    /// __str__
  std::string const str() const;


private:

  T val_;

};

//-----------------------------------------------------------------------------
template<typename T>
  Object::repr_t const Type<T>::repr() const
{
  std::stringstream ss;
  ss << val_;
  return repr_t(ss.str());
}

//-----------------------------------------------------------------------------
template<typename T>
  std::string const Type<T>::str() const
{
  std::stringstream ss;
  ss << val_;
  return ss.str();
}

} /* namespace icorne */
#endif /* __UFL_TYPE_H_ */
