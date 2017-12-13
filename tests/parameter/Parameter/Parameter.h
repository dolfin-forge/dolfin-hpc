#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/parameter/Parameter.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_Parameter )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_Parameter");
  {
    {
      Parameter P(bool(false));
      std::cout << P << "\n";
      Parameter Q(bool(true));
      std::cout << Q << "\n";
    }
    {
      Parameter P(int(0));
      std::cout << P << "\n";
      Parameter Q(int(1));
      std::cout << Q << "\n";
    }
    {
      Parameter P(uint(0u));
      std::cout << P << "\n";
      Parameter Q(uint(1u));
      std::cout << Q << "\n";
    }
    {
      Parameter P(real(0.0));
      std::cout << P << "\n";
      Parameter Q(real(1.0));
      std::cout << Q << "\n";
    }
    {
      Parameter P(std::string("0"));
      std::cout << P << "\n";
      Parameter Q(std::string("1"));
      std::cout << Q << "\n";
    }
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
