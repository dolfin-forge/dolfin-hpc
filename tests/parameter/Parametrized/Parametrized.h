#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/parameter/Parametrized.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_Parametrized )
  {
    Parametrized P0;
    ck_assert(P0.empty());
    P0.add("bool", bool());
    P0.add("int", int());
    P0.add("uint", uint());
    P0.add("real", real());
    P0.add("string", std::string());
    P0.disp();

    Parametrized P1;
    ck_assert(P1.empty());
    ck_assert(P0 != P1);
    P1.add("bool", bool());
    ck_assert(P0 != P1);
    P1.add("int", int());
    ck_assert(P0 != P1);
    P1.add("uint", uint());
    ck_assert(P0 != P1);
    P1.add("real", real());
    ck_assert(P0 != P1);
    P1.add("string", std::string());
    ck_assert(P0 == P1);
    P1.disp();

    Parametrized P2;
    ck_assert(P2.empty());
    ck_assert(P0 != P2);
    P2.add("bool", true);
    ck_assert(P0 != P2);
    P2.add("int", 1);
    ck_assert(P0 != P2);
    P2.add("uint", 1u);
    ck_assert(P0 != P2);
    P2.add("real", 1.0);
    ck_assert(P0 != P2);
    P2.add("string", "1");
    ck_assert(P0 != P2);
    P2.disp();

    P0 << P2;
    ck_assert(P0 == P2);
    ck_assert(P0 != P1);

    P1 >> P2;
    ck_assert(P1 == P2);
    ck_assert(P0 != P2);
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
