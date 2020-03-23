#include <dolfin_tests.h>
#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_CHECK

#include "File/File.h"
#include "STL/STL.h"
#include "VTK/VTK.h"
#include "Binary/Binary.h"

//-----------------------------------------------------------------------------
DOLFIN_SUITE_BEGIN(suite, "io")
{
  DOLFIN_TCASE_CREATE("File");
  DOLFIN_TCASE_ADD(test_File);

  DOLFIN_TCASE_CREATE("STL");
  DOLFIN_TCASE_ADD(test_STLMesh);

  DOLFIN_TCASE_CREATE("VTK");
  DOLFIN_TCASE_ADD(test_VTKMesh);

  DOLFIN_TCASE_CREATE("Binary");
  DOLFIN_TCASE_ADD(test_BinaryFile_Vector);
  DOLFIN_TCASE_ADD(test_BinaryFile_Mesh);
  DOLFIN_TCASE_ADD(test_BinaryFile_Poisson);
  DOLFIN_TCASE_ADD(test_BinaryFile_Function);
  DOLFIN_TCASE_ADD(test_BinaryFile_LabelList);
  DOLFIN_TCASE_ADD(test_BinaryFile_MeshFunction);
}
DOLFIN_SUITE_END
//-----------------------------------------------------------------------------
DOLFIN_CHECK_SUITE("dolfin/io", suite)
//-----------------------------------------------------------------------------

#endif
