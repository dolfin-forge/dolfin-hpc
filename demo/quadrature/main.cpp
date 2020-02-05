// Copyright (C) 2003-2005 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin.h>

#include <cstdlib>
#include <cstring>

using namespace dolfin;

int main(int argc, char** argv)
{
  if ( argc != 3 )
  {
    dolfin::cout << "Usage: dolfin-quadrature rule n' where rule is one of\n";
    dolfin::cout << "gauss, radau, lobatto, and n is the number of points\n";
    return 1;
  }

  int n = atoi(argv[2]);

  if ( strcmp(argv[1], "gauss") == 0 ) {

    GaussQuadrature q(n);
    q.disp();

  }
  else if ( strcmp(argv[1], "radau") == 0 ) {

    RadauQuadrature q(n);
    q.disp();

  }
  else if ( strcmp(argv[1], "lobatto") == 0 ) {

     LobattoQuadrature q(n);
     q.disp();

  }
  else {
    dolfin::cout << "Unknown quadrature rule.\n";
    return 1;
  }

  return 0;
}
