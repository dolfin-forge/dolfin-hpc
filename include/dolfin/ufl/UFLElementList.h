// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef  __UFL_ELEMENT_LIST_H
#define  __UFL_ELEMENT_LIST_H

#include <dolfin/ufl/UFLObject.h>

#include <map>
#include <string>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  ElementList
 *
 *  @brief  Provides the list of elements supported by UFL and specified in file
 *          ufl.elementlist of UFL version 2.1.1.
 *          Only cell types for simplicial meshes have been retained as
 *          dolfin-hpc does not support quadrilateral/hexahedral meshes.
 *
 */

class FiniteElementSpace;

class ElementList : protected std::map<Object::repr_t, FiniteElementSpace *>
{

  typedef std::pair<ufl::Object::repr_t, FiniteElementSpace *> ElementItem;

public:

  ///
  ElementList();

  //
  ElementList(ElementList const &other);

  ///
  ~ElementList();

  /// Add element in the list, nothing is done if it is already registered
  void add(Object::repr_t const& signature);

  /// Check if element is contained in the list
  bool has(Object::repr_t const& signature) const;

  /// Set the iterator to the begin of the list and return first element
  /// If the list is empty a null pointer is returned.
  FiniteElementSpace const * first() const;

  /// Set the iterator to the next element
  FiniteElementSpace const * next() const;

  /// Check if the end of the list is reached
  bool valid() const;

  /// Check the consistency of the list
  bool check() const;

  /// Display information on the list
  void disp() const;

private:

  mutable ElementList::const_iterator it_;

};

}

#endif /* __DOLFIN_UFL_ELEMENT_LIST_H */
