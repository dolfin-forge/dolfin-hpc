// Automatically generated code mapping element and dof map signatures
// to the corresponding ufc::finite_element and ufc::dofmap classes

#include <dolfin/log/log.h>

#include <cstring>

#include "ffc_Lagrange_1_1d.h"
#include "ffc_Lagrange_2_1d.h"
#include "ffc_Lagrange_1_2d.h"
#include "ffc_Lagrange_2_2d.h"
#include "ffc_Lagrange_1_3d.h"
#include "ffc_Lagrange_2_3d.h"
#include "ffc_Discontinuous_Lagrange_0_1d.h"
#include "ffc_Discontinuous_Lagrange_1_1d.h"
#include "ffc_Discontinuous_Lagrange_2_1d.h"
#include "ffc_Discontinuous_Lagrange_0_2d.h"
#include "ffc_Discontinuous_Lagrange_1_2d.h"
#include "ffc_Discontinuous_Lagrange_2_2d.h"
#include "ffc_Discontinuous_Lagrange_0_3d.h"
#include "ffc_Discontinuous_Lagrange_1_3d.h"
#include "ffc_Discontinuous_Lagrange_2_3d.h"
#include "ffc_Lagrange_1_2dVector.h"
#include "ffc_Lagrange_2_2dVector.h"
#include "ffc_Lagrange_1_3dVector.h"
#include "ffc_Lagrange_2_3dVector.h"
#include "ffc_Discontinuous_Lagrange_0_2dVector.h"
#include "ffc_Discontinuous_Lagrange_1_2dVector.h"
#include "ffc_Discontinuous_Lagrange_2_2dVector.h"
#include "ffc_Discontinuous_Lagrange_0_3dVector.h"
#include "ffc_Discontinuous_Lagrange_1_3dVector.h"
#include "ffc_Discontinuous_Lagrange_2_3dVector.h"
#include "ffc_Brezzi_Douglas_Marini_1_2d.h"

#include <dolfin/elements/ElementLibrary.h>

ufc::finite_element* dolfin::ElementLibrary::create_finite_element(const char* signature)
{
  if (strcmp(signature, "FiniteElement('Lagrange', interval, 1)") == 0)
    return new ffc_lagrange_1_1d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Lagrange', interval, 2)") == 0)
    return new ffc_lagrange_2_1d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Lagrange', triangle, 1)") == 0)
    return new ffc_lagrange_1_2d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Lagrange', triangle, 2)") == 0)
    return new ffc_lagrange_2_2d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Lagrange', tetrahedron, 1)") == 0)
    return new ffc_lagrange_1_3d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Lagrange', tetrahedron, 2)") == 0)
    return new ffc_lagrange_2_3d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', interval, 0)") == 0)
    return new ffc_discontinuous_lagrange_0_1d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', interval, 1)") == 0)
    return new ffc_discontinuous_lagrange_1_1d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', interval, 2)") == 0)
    return new ffc_discontinuous_lagrange_2_1d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', triangle, 0)") == 0)
    return new ffc_discontinuous_lagrange_0_2d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', triangle, 1)") == 0)
    return new ffc_discontinuous_lagrange_1_2d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', triangle, 2)") == 0)
    return new ffc_discontinuous_lagrange_2_2d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', tetrahedron, 0)") == 0)
    return new ffc_discontinuous_lagrange_0_3d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', tetrahedron, 1)") == 0)
    return new ffc_discontinuous_lagrange_1_3d_finite_element_0();
  if (strcmp(signature, "FiniteElement('Discontinuous Lagrange', tetrahedron, 2)") == 0)
    return new ffc_discontinuous_lagrange_2_3d_finite_element_0();
  if (strcmp(signature, "VectorElement(FiniteElement('Lagrange', triangle, 1), dim=2)") == 0)
    return new ffc_lagrange_1_2dvector_finite_element_1();
  if (strcmp(signature, "VectorElement(FiniteElement('Lagrange', triangle, 2), dim=2)") == 0)
    return new ffc_lagrange_2_2dvector_finite_element_1();
  if (strcmp(signature, "VectorElement(FiniteElement('Lagrange', tetrahedron, 1), dim=3)") == 0)
    return new ffc_lagrange_1_3dvector_finite_element_1();
  if (strcmp(signature, "VectorElement(FiniteElement('Lagrange', tetrahedron, 2), dim=3)") == 0)
    return new ffc_lagrange_2_3dvector_finite_element_1();
  if (strcmp(signature, "VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 0), dim=2)") == 0)
    return new ffc_discontinuous_lagrange_0_2dvector_finite_element_1();
  if (strcmp(signature, "VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 1), dim=2)") == 0)
    return new ffc_discontinuous_lagrange_1_2dvector_finite_element_1();
  if (strcmp(signature, "VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 2), dim=2)") == 0)
    return new ffc_discontinuous_lagrange_2_2dvector_finite_element_1();
  if (strcmp(signature, "VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 0), dim=3)") == 0)
    return new ffc_discontinuous_lagrange_0_3dvector_finite_element_1();
  if (strcmp(signature, "VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 1), dim=3)") == 0)
    return new ffc_discontinuous_lagrange_1_3dvector_finite_element_1();
  if (strcmp(signature, "VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 2), dim=3)") == 0)
    return new ffc_discontinuous_lagrange_2_3dvector_finite_element_1();
  if (strcmp(signature, "FiniteElement('Brezzi-Douglas-Marini', triangle, 1)") == 0)
    return new ffc_brezzi_douglas_marini_1_2d_finite_element_0();
  dolfin::error("Requested finite element '%s' has not been pregenerated.", signature);
  return 0;
}

ufc::dofmap* dolfin::ElementLibrary::create_dof_map(const char* signature)
{
  if (strcmp(signature, "FFC dofmap for FiniteElement('Lagrange', interval, 1)") == 0)
    return new ffc_lagrange_1_1d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Lagrange', interval, 2)") == 0)
    return new ffc_lagrange_2_1d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Lagrange', triangle, 1)") == 0)
    return new ffc_lagrange_1_2d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Lagrange', triangle, 2)") == 0)
    return new ffc_lagrange_2_2d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Lagrange', tetrahedron, 1)") == 0)
    return new ffc_lagrange_1_3d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Lagrange', tetrahedron, 2)") == 0)
    return new ffc_lagrange_2_3d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Discontinuous Lagrange', interval, 0)") == 0)
    return new ffc_discontinuous_lagrange_0_1d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Discontinuous Lagrange', interval, 1)") == 0)
    return new ffc_discontinuous_lagrange_1_1d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Discontinuous Lagrange', interval, 2)") == 0)
    return new ffc_discontinuous_lagrange_2_1d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Discontinuous Lagrange', triangle, 0)") == 0)
    return new ffc_discontinuous_lagrange_0_2d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Discontinuous Lagrange', triangle, 1)") == 0)
    return new ffc_discontinuous_lagrange_1_2d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Discontinuous Lagrange', triangle, 2)") == 0)
    return new ffc_discontinuous_lagrange_2_2d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Discontinuous Lagrange', tetrahedron, 0)") == 0)
    return new ffc_discontinuous_lagrange_0_3d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Discontinuous Lagrange', tetrahedron, 1)") == 0)
    return new ffc_discontinuous_lagrange_1_3d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Discontinuous Lagrange', tetrahedron, 2)") == 0)
    return new ffc_discontinuous_lagrange_2_3d_dofmap_0();
  if (strcmp(signature, "FFC dofmap for VectorElement(FiniteElement('Lagrange', triangle, 1), dim=2)") == 0)
    return new ffc_lagrange_1_2dvector_dofmap_1();
  if (strcmp(signature, "FFC dofmap for VectorElement(FiniteElement('Lagrange', triangle, 2), dim=2)") == 0)
    return new ffc_lagrange_2_2dvector_dofmap_1();
  if (strcmp(signature, "FFC dofmap for VectorElement(FiniteElement('Lagrange', tetrahedron, 1), dim=3)") == 0)
    return new ffc_lagrange_1_3dvector_dofmap_1();
  if (strcmp(signature, "FFC dofmap for VectorElement(FiniteElement('Lagrange', tetrahedron, 2), dim=3)") == 0)
    return new ffc_lagrange_2_3dvector_dofmap_1();
  if (strcmp(signature, "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 0), dim=2)") == 0)
    return new ffc_discontinuous_lagrange_0_2dvector_dofmap_1();
  if (strcmp(signature, "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 1), dim=2)") == 0)
    return new ffc_discontinuous_lagrange_1_2dvector_dofmap_1();
  if (strcmp(signature, "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 2), dim=2)") == 0)
    return new ffc_discontinuous_lagrange_2_2dvector_dofmap_1();
  if (strcmp(signature, "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 0), dim=3)") == 0)
    return new ffc_discontinuous_lagrange_0_3dvector_dofmap_1();
  if (strcmp(signature, "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 1), dim=3)") == 0)
    return new ffc_discontinuous_lagrange_1_3dvector_dofmap_1();
  if (strcmp(signature, "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 2), dim=3)") == 0)
    return new ffc_discontinuous_lagrange_2_3dvector_dofmap_1();
  if (strcmp(signature, "FFC dofmap for FiniteElement('Brezzi-Douglas-Marini', triangle, 1)") == 0)
    return new ffc_brezzi_douglas_marini_1_2d_dofmap_0();
  dolfin::error("Requested dofmap '%s' has not been pregenerated.", signature);
  return 0;
}

ufl::ElementList const dolfin::ElementLibrary::init_elements()
{
  ufl::ElementList list;
  list.add(ufl::Object::repr_t("FiniteElement('Lagrange', interval, 1)"));
  list.add(ufl::Object::repr_t("FiniteElement('Lagrange', interval, 2)"));
  list.add(ufl::Object::repr_t("FiniteElement('Lagrange', triangle, 1)"));
  list.add(ufl::Object::repr_t("FiniteElement('Lagrange', triangle, 2)"));
  list.add(ufl::Object::repr_t("FiniteElement('Lagrange', tetrahedron, 1)"));
  list.add(ufl::Object::repr_t("FiniteElement('Lagrange', tetrahedron, 2)"));
  list.add(ufl::Object::repr_t("FiniteElement('Discontinuous Lagrange', interval, 0)"));
  list.add(ufl::Object::repr_t("FiniteElement('Discontinuous Lagrange', interval, 1)"));
  list.add(ufl::Object::repr_t("FiniteElement('Discontinuous Lagrange', interval, 2)"));
  list.add(ufl::Object::repr_t("FiniteElement('Discontinuous Lagrange', triangle, 0)"));
  list.add(ufl::Object::repr_t("FiniteElement('Discontinuous Lagrange', triangle, 1)"));
  list.add(ufl::Object::repr_t("FiniteElement('Discontinuous Lagrange', triangle, 2)"));
  list.add(ufl::Object::repr_t("FiniteElement('Discontinuous Lagrange', tetrahedron, 0)"));
  list.add(ufl::Object::repr_t("FiniteElement('Discontinuous Lagrange', tetrahedron, 1)"));
  list.add(ufl::Object::repr_t("FiniteElement('Discontinuous Lagrange', tetrahedron, 2)"));
  list.add(ufl::Object::repr_t("VectorElement(FiniteElement('Lagrange', triangle, 1), dim=2)"));
  list.add(ufl::Object::repr_t("VectorElement(FiniteElement('Lagrange', triangle, 2), dim=2)"));
  list.add(ufl::Object::repr_t("VectorElement(FiniteElement('Lagrange', tetrahedron, 1), dim=3)"));
  list.add(ufl::Object::repr_t("VectorElement(FiniteElement('Lagrange', tetrahedron, 2), dim=3)"));
  list.add(ufl::Object::repr_t("VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 0), dim=2)"));
  list.add(ufl::Object::repr_t("VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 1), dim=2)"));
  list.add(ufl::Object::repr_t("VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 2), dim=2)"));
  list.add(ufl::Object::repr_t("VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 0), dim=3)"));
  list.add(ufl::Object::repr_t("VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 1), dim=3)"));
  list.add(ufl::Object::repr_t("VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 2), dim=3)"));
  list.add(ufl::Object::repr_t("FiniteElement('Brezzi-Douglas-Marini', triangle, 1)"));
  return list;
}

