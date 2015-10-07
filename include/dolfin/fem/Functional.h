// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-27
// Last changed: 2014-02-27

#ifndef __DOLFIN_FUNCTIONAL_H
#define __DOLFIN_FUNCTIONAL_H

#include <dolfin/fem/Form.h>

#include <ufc.h>

namespace dolfin
{

class Functional : public Form
{
public:

  /// Constructor
  Functional(Mesh& mesh);

  /// Destructor
  ~Functional();

private:

};

//--- INLINES -----------------------------------------------------------------

}

#endif
