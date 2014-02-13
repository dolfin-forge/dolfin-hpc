// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-02-06
// Last changed: 2014-02-06

#include <dolfin/function/SubFunction.h>

#include <dolfin/function/DiscreteFunction.h>
#include <dolfin/log/log.h>
#include <dolfin/log/LogStream.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
void SubFunction::disp() const
{
  cout << "SubFunction" << endl;
  cout << "-----------" << endl;

  // Begin indentation
  begin("");
  cout << "Index                 : " << this->index() << endl;
  skip();
  this->function().disp();
  // End indentation
  end();
  skip();
}

}
