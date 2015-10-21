#include <dolfin.h>

using namespace dolfin;

int main(int argc, char** argv)
{
  int ret = 0;

  //---------------------------------------------------------------------------
  Test t(argc, argv);

  {
    Mesh mesh0(t.args.mesh_file);
    std::string filename("save.xml");

    //-------------------------------------------------------------------------
    File file0("0"+filename);
    file0 << mesh0;

    Mesh mesh1("0"+filename);
    File file1("1"+filename);
    file1 << mesh1;

  }

  return ret;
}

