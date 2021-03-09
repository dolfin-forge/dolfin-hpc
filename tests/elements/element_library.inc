
#ifndef ELEMENT_LIBRARY_INC
#define ELEMENT_LIBRARY_INC

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

namespace ElementLibrary
{

constexpr char const * elements[] = {
  "FiniteElement('Lagrange', interval, 1)",
  "FiniteElement('Lagrange', interval, 2)",
  "FiniteElement('Lagrange', triangle, 1)",
  "FiniteElement('Lagrange', triangle, 2)",
  "FiniteElement('Lagrange', tetrahedron, 1)",
  "FiniteElement('Lagrange', tetrahedron, 2)",
  "FiniteElement('Discontinuous Lagrange', interval, 0)",
  "FiniteElement('Discontinuous Lagrange', interval, 1)",
  "FiniteElement('Discontinuous Lagrange', interval, 2)",
  "FiniteElement('Discontinuous Lagrange', triangle, 0)",
  "FiniteElement('Discontinuous Lagrange', triangle, 1)",
  "FiniteElement('Discontinuous Lagrange', triangle, 2)",
  "FiniteElement('Discontinuous Lagrange', tetrahedron, 0)",
  "FiniteElement('Discontinuous Lagrange', tetrahedron, 1)",
  "FiniteElement('Discontinuous Lagrange', tetrahedron, 2)",
  "VectorElement(FiniteElement('Lagrange', triangle, 1), dim=2)",
  "VectorElement(FiniteElement('Lagrange', triangle, 2), dim=2)",
  "VectorElement(FiniteElement('Lagrange', tetrahedron, 1), dim=3)",
  "VectorElement(FiniteElement('Lagrange', tetrahedron, 2), dim=3)",
  "VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 0), dim=2)",
  "VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 1), dim=2)",
  "VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 2), dim=2)",
  "VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 0), dim=3)",
  "VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 1), dim=3)",
  "VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 2), dim=3)",
  "FiniteElement('Brezzi-Douglas-Marini', triangle, 1)"
};

constexpr char const * dofmaps[] = {
  "FFC dofmap for FiniteElement('Lagrange', interval, 1)",
  "FFC dofmap for FiniteElement('Lagrange', interval, 2)",
  "FFC dofmap for FiniteElement('Lagrange', triangle, 1)",
  "FFC dofmap for FiniteElement('Lagrange', triangle, 2)",
  "FFC dofmap for FiniteElement('Lagrange', tetrahedron, 1)",
  "FFC dofmap for FiniteElement('Lagrange', tetrahedron, 2)",
  "FFC dofmap for FiniteElement('Discontinuous Lagrange', interval, 0)",
  "FFC dofmap for FiniteElement('Discontinuous Lagrange', interval, 1)",
  "FFC dofmap for FiniteElement('Discontinuous Lagrange', interval, 2)",
  "FFC dofmap for FiniteElement('Discontinuous Lagrange', triangle, 0)",
  "FFC dofmap for FiniteElement('Discontinuous Lagrange', triangle, 1)",
  "FFC dofmap for FiniteElement('Discontinuous Lagrange', triangle, 2)",
  "FFC dofmap for FiniteElement('Discontinuous Lagrange', tetrahedron, 0)",
  "FFC dofmap for FiniteElement('Discontinuous Lagrange', tetrahedron, 1)",
  "FFC dofmap for FiniteElement('Discontinuous Lagrange', tetrahedron, 2)",
  "FFC dofmap for VectorElement(FiniteElement('Lagrange', triangle, 1), dim=2)",
  "FFC dofmap for VectorElement(FiniteElement('Lagrange', triangle, 2), dim=2)",
  "FFC dofmap for VectorElement(FiniteElement('Lagrange', tetrahedron, 1), dim=3)",
  "FFC dofmap for VectorElement(FiniteElement('Lagrange', tetrahedron, 2), dim=3)",
  "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 0), dim=2)",
  "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 1), dim=2)",
  "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', triangle, 2), dim=2)",
  "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 0), dim=3)",
  "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 1), dim=3)",
  "FFC dofmap for VectorElement(FiniteElement('Discontinuous Lagrange', tetrahedron, 2), dim=3)",
  "FFC dofmap for FiniteElement('Brezzi-Douglas-Marini', triangle, 1)"
};

inline ufc::finite_element * create_finite_element(const char* signature)
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

inline ufc::dofmap* create_dof_map(const char* signature)
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

} // namespace ElementLibrary

#endif
