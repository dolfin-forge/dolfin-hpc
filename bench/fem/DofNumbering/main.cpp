#include <dolfin.h>

#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/fem/DofNumbering.h>
#include <dolfin/fem/FiniteElementSpace.h>

#include <dolfin/mesh/IntervalCell.h>
#include <dolfin/mesh/TriangleCell.h>
#include <dolfin/mesh/TetrahedronCell.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
std::string dofmap_of(ufl::FiniteElementSpace& e)
{
  return DofMap::make_signature(e.repr());
}
//-----------------------------------------------------------------------------
void test(DofNumbering& dofnum, UFCMesh& ufc_mesh, ufc::dofmap& ufc_dofmap)
{
  dofnum.init(ufc_mesh, ufc_dofmap);
  dofnum.build(ufc_mesh, ufc_dofmap);
  dofnum.disp();
  uint * mapping = NULL;
  uint mapping_size = 0;
  dofnum.pretabulate(mapping, mapping_size);
  dofnum.cache();
  dolfin_assert(dofnum.map_size() == mapping_size);
  Mesh& mesh = const_cast<Mesh&>(*ufc_mesh.mesh);
  DofMap dm(mesh, ufc_dofmap, false);
  uint * dofs = new uint[ufc_dofmap.local_dimension()];
  CellIterator c(mesh);
  UFCCell ufc_cell(*c);
  uint ii = 0;
  for (; !c.end(); ++c)
  {
    ufc_cell.update(*c);
    dm.tabulate_dofs(dofs, ufc_cell, *c);
    for (uint i = 0; i < ufc_dofmap.local_dimension(); ++i, ++ii)
    {
      if(mapping[ii] != dofs[i])
      {
        error("Invalid dofs mapping");
      }
    }
  }
  delete[] dofs;
  delete[] mapping;
}

//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  dolfin_init(argc, argv);
  //---------------------------------------------------------------------------
  {

    // 1D: broken in parallel
		if(MPI::size() == 1)
    {
      begin("1D");
      Mesh mesh("../../data/meshes/intervalm.bin");
      mesh.distdata().check(true);
      UFCMesh ufc_mesh(mesh);
      ufc::dofmap * ufc_dofmap;
      uint const tdim = mesh.topology().dim();
      //-----------------------------------------------------------------------
      {
        message("DG0s");
        ufl::FiniteElement DG0s(ufl::Family::DG, mesh.type(), 0);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(DG0s));
        DG0sNumbering dnDG0s;
        test(dnDG0s, ufc_mesh, *ufc_dofmap);;
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("CG1s");
        ufl::FiniteElement CG1s(ufl::Family::CG, mesh.type(), 1);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(CG1s));
        CG1sNumbering dnCG1s;
        test(dnCG1s, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("CG2s");
        ufl::FiniteElement CG2s(ufl::Family::CG, mesh.type(), 2);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(CG2s));
        Parallel0Numbering dnCG2s;
        test(dnCG2s, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      end();
      skip();
    }

    // 2D
    {
      begin("2D");
      Mesh mesh("../../data/meshes/trianglem.bin");
      mesh.distdata().check(true);
      UFCMesh ufc_mesh(mesh);
      ufc::dofmap * ufc_dofmap;
      uint const tdim = mesh.topology().dim();
      //-----------------------------------------------------------------------
      {
        message("DG0s");
        ufl::FiniteElement DG0s(ufl::Family::DG, mesh.type(), 0);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(DG0s));
        DG0sNumbering dnDG0s;
        test(dnDG0s, ufc_mesh, *ufc_dofmap);;
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("DG0v");
        ufl::VectorElement DG0v(ufl::Family::DG, mesh.type(), 0, tdim);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(DG0v));
        DG0vNumbering dnDG0v;
        test(dnDG0v, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("CG1s");
        ufl::FiniteElement CG1s(ufl::Family::CG, mesh.type(), 1);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(CG1s));
        CG1sNumbering dnCG1s;
        test(dnCG1s, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("CG1v");
        ufl::VectorElement CG1v(ufl::Family::CG, mesh.type(), 1, tdim);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(CG1v));
        CG1vNumbering dnCG1v;
        test(dnCG1v, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("CG2s");
        ufl::FiniteElement CG2s(ufl::Family::CG, mesh.type(), 2);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(CG2s));
        Parallel0Numbering dnCG2s;
        test(dnCG2s, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("CG2v");
        ufl::VectorElement CG2v(ufl::Family::CG, mesh.type(), 2, tdim);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(CG2v));
        Parallel1Numbering dnCG2v;
        test(dnCG2v, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      end();
      skip();
    }

    // 3D
    {
      begin("3D");
      Mesh mesh("../../data/meshes/tetrahedronm.bin");
      mesh.distdata().check(true);
      UFCMesh ufc_mesh(mesh);
      ufc::dofmap * ufc_dofmap;
      uint const tdim = mesh.topology().dim();
      //-----------------------------------------------------------------------
      {
        message("DG0s");
        ufl::FiniteElement DG0s(ufl::Family::DG, mesh.type(), 0);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(DG0s));
        DG0sNumbering dnDG0s;
        test(dnDG0s, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("DG0v");
        ufl::VectorElement DG0v(ufl::Family::DG, mesh.type(), 0, tdim);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(DG0v));
        DG0vNumbering dnDG0v;
        test(dnDG0v, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("CG1s");
        ufl::FiniteElement CG1s(ufl::Family::CG, mesh.type(), 1);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(CG1s));
        CG1sNumbering dnCG1s;
        test(dnCG1s, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("CG1v");
        ufl::VectorElement CG1v(ufl::Family::CG, mesh.type(), 1, tdim);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(CG1v));
        CG1vNumbering dnCG1v;
        test(dnCG1v, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("CG2s");
        ufl::FiniteElement CG2s(ufl::Family::CG, mesh.type(), 2);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(CG2s));
        Parallel0Numbering dnCG2s;
        test(dnCG2s, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      {
        message("CG2v");
        ufl::VectorElement CG2v(ufl::Family::CG, mesh.type(), 2, tdim);
        ufc_dofmap = ElementLibrary::create_dof_map(dofmap_of(CG2v));
        Parallel1Numbering dnCG2v;
        test(dnCG2v, ufc_mesh, *ufc_dofmap);
        delete ufc_dofmap;
      }
      //-----------------------------------------------------------------------
      end();
      skip();
    }

  }
  //---------------------------------------------------------------------------
  dolfin_finalize();
  return 0;
}

