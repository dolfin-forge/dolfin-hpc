
#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include <check.h>
#include <dolfin/ufl/UFLFamily.h>


using ufl::Family;

//-----------------------------------------------------------------------------
START_TEST( test_UFL_Family )
{
  int init_failed = 0;

  std::vector<Family::Type> v;
  v.push_back(Family::ARG);
  v.push_back(Family::AW);
  v.push_back(Family::BDFM);
  v.push_back(Family::BDM);
  v.push_back(Family::CR);
  v.push_back(Family::DG);
  v.push_back(Family::HER);
  v.push_back(Family::CG);
  v.push_back(Family::MTW);
  v.push_back(Family::MOR);
  v.push_back(Family::N1curl);
  v.push_back(Family::N2curl);
  v.push_back(Family::RT);
  v.push_back(Family::BQ);
  v.push_back(Family::B);
  v.push_back(Family::Q);
  v.push_back(Family::R);
  v.push_back(Family::U);


  for (std::vector<Family::Type>::const_iterator it = v.begin();
       it != v.end(); ++it)
  {
    Family f(*it);
    f.display();
  }
  std::cout << std::endl;

  fail_unless( init_failed == 0 );


  fail_unless( init_failed == 0 );
}END_TEST
//-----------------------------------------------------------------------------


#endif
