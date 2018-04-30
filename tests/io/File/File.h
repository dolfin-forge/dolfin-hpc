#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/io/File.h>
#include <dolfin/parameter/parameters.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_File )
  {
      std::string basename("basename");
      dolfin_set("output_format", "binary");
      message("filename::binary = %s", File::filename(basename).c_str());
      ck_assert(File::filename(basename) == (basename + ".bin"));
      dolfin_set("output_format", "vtk");
      message("filename::vtk    = %s", File::filename(basename).c_str());
      ck_assert(File::filename(basename) == (basename + ".pvd"));
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
