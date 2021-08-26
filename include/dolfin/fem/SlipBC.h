// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.

// Existing code for Dirichlet BC is used

#ifndef __DOLFIN_SLIPBC_H
#define __DOLFIN_SLIPBC_H

#include <dolfin/fem/BoundaryCondition.h>
#include <dolfin/fem/NodeNormal.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>

namespace dolfin
{

class DofMap;
class Form;
class Function;
class Mesh;
class ScratchSpace;
class SubDomain;

class SlipBC : public BoundaryCondition
{
public:
  /// Create boundary condition for sub domain
  SlipBC( Mesh & mesh, SubDomain const & sub_domain );

  /// Create boundary condition for sub domain given normals
  SlipBC( Mesh & mesh, SubDomain const & sub_domain, NodeNormal & normals );

  /// Create sub system boundary condition for sub domain
  SlipBC( Mesh &            mesh,
          SubDomain const & sub_domain,
          SubSystem const & sub_system );

  /// Destructor
  ~SlipBC() override;

  /// Apply boundary condition to linear system
  void apply( GenericMatrix &      A,
              GenericVector &      b,
              BilinearForm const & form ) override;

  /// Apply boundary condition to non linear system
  void apply( GenericMatrix &       A,
              GenericVector &       b,
              GenericVector const & x,
              BilinearForm const &  form ) override;

  /// Apply boundary condition to off-diagonal block in linear system
  void zero( GenericMatrix &      A,
              BilinearForm const & form );

  auto normal() -> BoundaryNormal &
  {
    return *node_normal;
  }

private:
  inline void sync( Time const & ) override
  { /* No-op */
  }

  void applySlipBC_P1( GenericMatrix &      A,
                       GenericVector &      b,
                       BilinearForm const & form,
                       ScratchSpace &       scratch );

  void applySlipBC( GenericMatrix &      A,
                    GenericVector &      b,
                    BilinearForm const & form,
                    ScratchSpace &       scratch );

  void applyNodeBC( GenericMatrix &               A,
                    GenericVector &               b,
                    Mesh const &                  mesh,
                    size_t const                  node,
                    std::vector< size_t > const & udofs,
                    std::vector< size_t > const & ndofs );

  void applyZeroSlipBC( GenericMatrix &      A,
                    BilinearForm const & form,
                    ScratchSpace &       scratch );

  void applyZeroNodeBC( GenericMatrix &               A,
                    Mesh const &                  mesh,
                    size_t const                  node,
                    std::vector< size_t > const & udofs,
                    std::vector< size_t > const & ndofs );

  // Initialize sub domain markers
  void init( SubDomain const & sub_domain );

  // The mesh
  Mesh & mesh;

  // Node normal and tangents
  NodeNormal * node_normal;
  bool         node_normal_local;

  Matrix * As;
  bool     As_local;

  // Local data structures for assembly
  _ordered_set< size_t > row_indices;
  std::vector< real >    a[3];          // local lhs extracted from A
  std::vector< real >    a_slip_row[3]; // local lhs row after slip enforcement
  std::vector< size_t >  a_col_indices[3]; // non-zero indices per row
  real                   l[3];             // local rhs extracted from b
  real                   l_slip[3];        // local rhs after slip enforcement
  real                   basis_[3][3];     // local basis (n, tau1, tau2)
  size_t                 max[3];           // maximum component
  size_t                 row[3];           // row reordering
};

} /* namespace dolfin */

#endif /* __DOLFIN_SLIPBC_H */
