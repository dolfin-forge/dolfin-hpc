#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/parameter/Parameter.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Parameter )
  {
    {
      Parameter P(bool(false));
      std::cout << P << "\n";
      Parameter Q(bool(true));
      std::cout << Q << "\n";
    }
    //---
    {
      Parameter P(int(0));
      std::cout << P << "\n";
      Parameter Q(int(1));
      std::cout << Q << "\n";
    }
    //---
    {
      Parameter P(uint(0u));
      std::cout << P << "\n";
      Parameter Q(uint(1u));
      std::cout << Q << "\n";
    }
    //---
    {
      Parameter P(real(0.0));
      std::cout << P << "\n";
      Parameter Q(real(1.0));
      std::cout << Q << "\n";
    }
    //---
    {
      Parameter P(std::string("0"));
      std::cout << P << "\n";
      Parameter Q(std::string("1"));
      std::cout << Q << "\n";
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
