#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/io/BinaryFile.h>

#include <dolfin/la/Vector.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>

#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/mesh/UnitInterval.h>
#include <dolfin/ufl/UFLFiniteElement.h>

#include <dolfin/dolfin.h>
#include "../../demo/pde/poisson/Poisson.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_BinaryFile_Vector )
  {
    //--------------------------------------------------------------------------
    // (Generic)Vector
    {
      Vector v1(50);

      {
        v1.zero();

        real *data = new real[50];

        for (int i = 0; i < 50; i++)
        {
          data[i] = 3.1415;
        }

        v1.set(data);
        delete[] data;
      }

      {
        BinaryFile f1("Vector.bin");
        f1 << v1;
      }

      Vector v2;

      {
        BinaryFile f2("Vector.bin");
        f2 >> v2;
      }

      ck_assert( v1.size() == v2.size() );
      ck_assert( v1[0] == v2[0] );
      ck_assert( v1[v1.size()/2] == v2[v2.size()/2] );

    }
  }
DOLFIN_END_TEST
    //--------------------------------------------------------------------------
    // Mesh
DOLFIN_START_TEST( test_BinaryFile_Mesh )
  {
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
        BinaryFile f1("mesh.bin");
        f1 << mesh1;
      }

      Mesh mesh2(*type, space);
      {
        BinaryFile f2("mesh.bin");
        f2 >> mesh2;
      }

      ck_assert( mesh1 == mesh2 );
    }
  }
DOLFIN_END_TEST
    //--------------------------------------------------------------------------
    // example from demo/pde/poisson
DOLFIN_START_TEST( test_BinaryFile_Poisson )
  {
    {
      //------------------------------------------------------------------------

      // Source term
      struct Source : public Value<Source>
      {
        void eval(real * value, const real* x) const
        {
          real dx = x[0] - 0.5;
          real dy = x[1] - 0.5;
          value[0] = 500.0*std::exp(-(dx*dx + dy*dy)/0.02);
        }
      };

      // Neumann boundary condition
      struct Flux : public Value<Flux>
      {
        void eval(real * value, const real* x) const
        {
          if (x[0] > DOLFIN_EPS)
            value[0] = 25.0*std::sin(5.0*DOLFIN_PI*x[1]);
          else
            value[0] = 0.0;
        }
      };

      //------------------------------------------------------------------------

      // Sub domain for Dirichlet boundary condition
      struct DirichletBoundary : public SubDomain
      {
        bool inside(const real* x, bool on_boundary) const
        {
          return x[0] < DOLFIN_EPS && on_boundary;
        }
      };

      //------------------------------------------------------------------------

      // Create mesh
      // Mesh mesh("./demo/pde/poisson/UnitSquareMesh_32x32.xml");
      UnitSquare mesh(32,32);

      // Create coefficients
      Analytic<Source>  f(mesh);
      Analytic<Flux>    g(mesh);

      // Create boundary condition
      Constant u0(0.0);
      DirichletBoundary boundary;
      DirichletBC bc(u0, mesh, boundary);

      // Define PDE
      Poisson::BilinearForm a(mesh);
      Poisson::LinearForm   L(mesh, f, g);

      // Solve PDE
      Matrix A;
      Vector b;

      a.assemble(A, true);
      L.assemble(b, true);
      bc.apply(A, b, a);


      Function u(a.trial_space());
      KrylovSolver solver(bicgstab, bjacobi);

      solver.solve(A, u.vector(), b);
      u.sync();

      u.disp();

      {
        BinaryFile fo("u.bin");
        fo << u;
      }

      message("vector l2  norm: %e", u.vector().norm());
      message("vector inf norm: %e", u.vector().max());

      Function v(a.trial_space());

      {
        BinaryFile fi("u000000.bin");
        fi >> v;
      }
    }
  }
DOLFIN_END_TEST
    //--------------------------------------------------------------------------
    // Function
DOLFIN_START_TEST( test_BinaryFile_Function )
  {
    {
      // taken from the function test case
      // UnitInterval m0(16);
      // ufl::FiniteElement DG0(ufl::Family::DG, m0.type(), 0);
      // FiniteElementSpace Vh0(m0, DG0);
      // Function U0(Vh0);

      // Mesh mesh("./demo/pde/poisson/UnitSquareMesh_32x32.xml");
      ufl::Family f(ufl::Family::CG);
      ufl::Domain dom(ufl::Domain::triangle);
      ufl::Cell cell(dom);
      // Mesh mesh;
      // CellType * ctype = CellType::create(cell);
      // ctype->create_reference_cell(mesh);
      // UnitCube mesh(32,32,32);
      UnitDisk mesh(32);
      ufl::FiniteElement uflfem(ufl::Family::CG, cell, 1);
      dolfin::FiniteElement fem(uflfem);
      ufc::dofmap * dmap = ElementLibrary::create_dof_map( DofMap::make_signature(fem.signature()));
      DofMap dm(mesh, *dmap, true);

      FiniteElementSpace femspace(mesh, fem, dm, false);
      Function U0(femspace);
      U0.vector() = 1.;

      U0.disp();

      {
        BinaryFile f1("function.bin");
        f1 << U0;
      }

      Function U1(femspace);
      U1.vector() = 0.;

      {
        BinaryFile f2("function000000.bin");
        f2 >> U1;
      }

      // U0.disp();
      // U1.disp();
      // ck_assert( U0 == U1 );
    }
  }
DOLFIN_END_TEST
    //--------------------------------------------------------------------------
    // LabelList<Function>
DOLFIN_START_TEST( test_BinaryFile_LabelList )
  {
    {
      // taken from the function test case
      UnitInterval m0(16);
      ufl::FiniteElement DG0(ufl::Family::DG, m0.type(), 0);
      FiniteElementSpace Vh0(m0, DG0);
      Function U0(Vh0);

      U0.disp();

      LabelList<Function> tmp1(1, Label<Function>(U0, "U0"));
      BinaryFile f1("LabelList.bin");
      f1 << tmp1;
      Function U1(Vh0);
      LabelList<Function> tmp2(1, Label<Function>(U1, "U0"));
      BinaryFile f2("LabelList000000.bin");
      f2 >> tmp2;
      // ck_assert( tmp1 == tmp2 );
    }
  }
DOLFIN_END_TEST
    //--------------------------------------------------------------------------
    // MeshFunction<real>
DOLFIN_START_TEST( test_BinaryFile_MeshFunction )
  {
    {
      // TODO
      // ck_assert( 3 == (2+1));
    }

  }
DOLFIN_END_TEST
//------------------------------------------------------------------------------

#endif
