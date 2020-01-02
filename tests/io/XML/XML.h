#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/io/XMLFile.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_XMLFile )
  {
    //-------------------------------------------------------------------------
    {
      CellType * type = CellType::create(CellType::triangle);
      EuclideanSpace space(2);
      Mesh mesh1(*type, space);

      {
        real coords[4][2] = {{0.0,0.0},
                             {1.0,0.0},
                             {1.0,1.0},
                             {0.0,1.0}};

        uint cells[2][3] = {{0, 1, 2},
                            {0, 2, 3}};

        MeshEditor editor(mesh1, CellType::triangle, 2 );
        editor.init_vertices(4);
        editor.init_cells(2);
        editor.add_vertex(0, coords[0]);
        editor.add_vertex(1, coords[1]);
        editor.add_vertex(2, coords[2]);
        editor.add_vertex(3, coords[3]);
        editor.add_cell(0, cells[0]);
        editor.add_cell(1, cells[1]);
        editor.close();
      }

      {
        XMLFile f1("mesh.xml");
        f1 << mesh1;
      }

      Mesh mesh2(*type, space);
      {
        XMLFile f2("mesh.xml");
        f2 >> mesh2;
      }
      ck_assert( mesh1 == mesh2 );
    }
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_XMLMesh )
  {
  }
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
