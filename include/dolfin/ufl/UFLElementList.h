// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef  __UFL_ELEMENT_LIST_H_
#define  __UFL_ELEMENT_LIST_H_

#include <dolfin/ufl/UFLDomain.h>
#include <dolfin/ufl/UFLFamily.h>

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>

#include <iomanip>
#include <map>
#include <set>
#include <sstream>

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

class ElementList
{

public:

  /// Meyers singleton
  static ElementList const& Supported()
  {
    static ElementList instance_;
    return instance_;
  }

private:

  ///
  ElementList();

  ///
  ~ElementList();

};

}

#endif /* __UFL_ELEMENT_LIST_H_ */
