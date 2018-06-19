#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/io/File.h>
#include <dolfin/parameter/parameters.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_File )
  {
    std::string basename("basename");
    std::string filename;
    dolfin_set("output_format", "binary");
		filename = File::filename(basename);
    message("filename::binary = %s", filename.c_str());
    ck_assert(filename == (basename + ".bin"));
    filename = File::filename(basename, "binary");
    message("filename::binary = %s", filename.c_str());
    dolfin_set("output_format", "vtk");
    filename = File::filename(basename);
    message("filename::vtk    = %s", filename.c_str());
    ck_assert(filename == (basename + ".pvd"));
    filename = File::filename(basename, "vtk");
    message("filename::vtk    = %s", filename.c_str());
    ck_assert(filename == (basename + ".pvd"));
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
