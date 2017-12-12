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
      Parameter P(bool());
    }
    {
      Parameter P(int());
    }
    {
      Parameter P(uint());
    }
    {
      Parameter P(real());
    }
    {
      Parameter P(std::string());
    }
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
