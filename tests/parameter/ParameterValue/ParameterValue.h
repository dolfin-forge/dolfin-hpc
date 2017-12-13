#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/parameter/ParameterValue.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
START_TEST( test_ParameterValue )
{
  int init_failed = 0;
  Test T;
  //---
  T.begin("test_ParameterValue");
  {
    message("bool");
    {
      parameter<bool> P;
      message("P = %u", bool(P));
      ck_assert(static_cast<bool>(P) == bool());
      parameter<bool> Q(false);
      message("Q = %u", bool(P));
      ck_assert(static_cast<bool>(Q) == false);
      ck_assert(P == Q);
      parameter<bool> R(true);
      message("R = %u", bool(P));
      ck_assert(static_cast<bool>(R) == true);
      ck_assert(P != R);
      ck_assert(Q != R);
      parameter<bool> S(R);
      ck_assert(R == S);
    }
    //---
    message("uint");
    {
      parameter<uint> P;
      message("P = %u", uint(P));
      ck_assert(static_cast<uint>(P) == uint());
      parameter<uint> Q(0);
      message("Q = %u", uint(P));
      ck_assert(static_cast<uint>(Q) == 0);
      ck_assert(P == Q);
      parameter<uint> R(1);
      message("E = %u", uint(P));
      ck_assert(static_cast<uint>(R) == 1);
      ck_assert(P != R);
      ck_assert(Q != R);
      static_cast<uint&>(Q) += 1;
      parameter<uint> S(Q);
      message("S = %u", uint(S));
      ck_assert(R == S);
    }
    //---
    message("int");
    {
      parameter<int> P;
      message("P = %+d", int(P));
      ck_assert(static_cast<int>(P) == int());
      parameter<int> Q(0);
      message("Q = %+d", int(Q));
      ck_assert(static_cast<int>(Q) == 0);
      ck_assert(P == Q);
      parameter<int> R(1);
      message("R = %+d", int(R));
      ck_assert(static_cast<int>(R) == 1);
      ck_assert(P != R);
      ck_assert(Q != R);
      parameter<int> S(-1);
      message("S = %+d", int(S));
      ck_assert(static_cast<int>(S) == -1);
      ck_assert(P != S);
      ck_assert(Q != S);
      ck_assert(R != S);
      static_cast<int&>(S) += 2;
      parameter<int> T(S);
      message("T = %+d", int(T));
      ck_assert(R == S);
      ck_assert(R == T);
    }
    //---
    message("real");
    {
      parameter<real> P;
      message("P = %+e", real(P));
      ck_assert(static_cast<real>(P) == real());
      parameter<real> Q(0.);
      message("Q = %+e", real(Q));
      ck_assert(static_cast<real>(Q) == 0.);
      ck_assert(P == Q);
      parameter<real> R(1);
      message("R = %+e", real(R));
      ck_assert(static_cast<real>(R) == 1.);
      ck_assert(P != R);
      ck_assert(Q != R);
      parameter<real> S(-1.);
      message("S = %+e", real(S));
      ck_assert(static_cast<real>(S) == -1.);
      ck_assert(P != S);
      ck_assert(Q != S);
      ck_assert(R != S);
      static_cast<real&>(S) = 1.;
      parameter<real> T(S);
      message("T = %+e", real(T));
      ck_assert(R == S);
      ck_assert(R == T);
    }
    //---
    message("string");
    {
      parameter<std::string> P;
      message("P = %s", std::string(P).c_str());
      ck_assert(static_cast<std::string>(P) == std::string());
      parameter<std::string> Q("");
      message("Q = %s", std::string(Q).c_str());
      ck_assert(static_cast<std::string>(Q) == "");
      ck_assert(P == Q);
      parameter<std::string> R("1");
      message("R = %s", std::string(P).c_str());
      ck_assert(static_cast<std::string>(R) == "1");
      ck_assert(P != R);
      ck_assert(Q != R);
      parameter<std::string> S("11");
      message("S = %s", std::string(P).c_str());
      ck_assert(static_cast<std::string>(S) == "11");
      ck_assert(P != S);
      ck_assert(Q != S);
      ck_assert(R != S);
      static_cast<std::string&>(R) += "1";
      parameter<std::string> T(S);
      message("T = %s", std::string(P).c_str());
      ck_assert(R == S);
      ck_assert(R == T);
    }
  }
  T.end();
  //---
  ck_assert( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------

#endif
