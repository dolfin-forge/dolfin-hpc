// Copyright (C) 2013 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  
// Last changed:

#ifndef __UFC_CELL_INTEGRAL_H
#define __UFC_CELL_INTEGRAL_H

#include <iomanip>
#include <dolfin/common/types.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/function/Function.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/quadrature/UFCReferenceCell.h>
#include <dolfin/quadrature/QuadratureRule.h>

namespace dolfin
{

  class UFCCellIntegral
  {
  public:
    
    /// Constructor
    UFCCellIntegral();

    /// Destructor
    virtual ~UFCCellIntegral();

    ///Tabulate the tensor for the contribution from a local cell
    virtual void tabulate_tensor(const UFC& ufc, Cell& cell,
        const Array<Function*>& coefficients,
        const QuadratureRule& q, unsigned int** dofs) const;
  };

  inline UFCCellIntegral::UFCCellIntegral()
  {}

  inline UFCCellIntegral::~UFCCellIntegral()
  {}

  inline void UFCCellIntegral::tabulate_tensor(const UFC& ufc, Cell& cell, 
      const Array<Function*>& coefficients,
      const QuadratureRule& q,
      unsigned int** dofs) const
  {
    // Extract vertex coordinates
    const double * const * x = ufc.cell.coordinates;
    
    UFCReferenceCell ref_cell(cell);
    const double * const * x_ref = ref_cell.coordinates;

    const unsigned int dim = ufc.mesh.topological_dimension;
    const unsigned int n_vertices = (dim == 1 ? 2 : (dim == 2 ? 3 : 4));
    std::vector<std::vector<double> > J (2, std::vector<double>(2,0.));
    J[0][0] = x[1][0] - x[0][0];
    J[0][1] = x[2][0] - x[0][0];
    J[1][0] = x[1][1] - x[0][1];
    J[1][1] = x[2][1] - x[0][1];
    // Compute determinant of Jacobian
    double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
      
    // Compute inverse of Jacobian
    std::vector<std::vector<double> > Jinv (2, std::vector<double>(2,0.));
    Jinv[0][0] =  J[1][1] / detJ;
    Jinv[0][1] = -J[0][1] / detJ;
    Jinv[1][0] = -J[1][0] / detJ;
    Jinv[1][1] =  J[0][0] / detJ;
    
    // Set scale factor
    const double det = std::abs(detJ);
    
    unsigned int tensor_rank = ufc.form.rank();
    std::vector<unsigned int> n_test(tensor_rank);
    for(unsigned int i=0; i<n_test.size(); ++i)
      n_test[i] = ufc.local_dimensions[i];
    
    // Reset values of the element tensor block
    for (unsigned int j = 0; j < n_test[0]; j++)
    {
      if(tensor_rank == 2)
      {
        for (unsigned int k = 0; k < n_test[0]; k++)
        {
          ufc.A[j*n_test[0] + k] = 0;
        }// end loop over 'k'
      }
      else if(tensor_rank == 1)
        ufc.A[j] = 0;
      else
        std::cout << "not implemented" << std::endl;
        
    }// end loop over 'j'
    
    const std::vector<real>& weights = q.get_weights();
    const std::vector<real*>& q_points = q.get_points();
    const std::vector<real*>& ref_points = q.get_reference_points();

    std::vector<real*> real_points (q_points.size());
    for(unsigned int i = 0; i<real_points.size(); ++i)
    {
        real_points[i] = new real[dim];
        for(unsigned int d=0; d<dim; ++d)
        {
          real_points[i][d] = 0.;
          for(unsigned int e=0; e<dim; ++e)
            real_points[i][d] += J[d][e] * q_points[i][e];
          real_points[i][d] += x[0][d];
        }
    }

    std::vector<std::vector<real*> > q_coefficients (coefficients.size());
    for(unsigned int i = 0; i<q_coefficients.size(); ++i)
    {
      q_coefficients[i].resize(q.size());
      for(unsigned int q = 0; q<q_coefficients[i].size(); ++q)
      {
        q_coefficients[i][q] = new real[1];
        coefficients[i]->eval(q_coefficients[i][q], real_points[q]);
      }
    }

    std::vector<std::vector<std::vector<real*> > > phi_values(tensor_rank);
    std::vector<std::vector<std::vector<real*> > > phi_grads(tensor_rank);

    std::vector<std::vector<std::vector<real*> > > ref_phi_values(tensor_rank);
    std::vector<std::vector<std::vector<real*> > > ref_phi_grads(tensor_rank);

    for (unsigned int i = 0; i<tensor_rank; ++i)
    {
      phi_values[i].resize(n_test[i]);
      phi_grads[i].resize(n_test[i]);
      for (unsigned int j = 0; j<n_test[i]; ++j)
      {
        phi_values[i][j].resize(q.size());
        phi_grads[i][j].resize(q.size());
        for (unsigned int qp = 0; qp<q.size(); ++qp)
        {
          phi_values[i][j][qp] = new real[1];
          phi_grads[i][j][qp] = new real[dim];
          ufc.finite_elements[i]->evaluate_basis(j, phi_values[i][j][qp], q_points[qp], ref_cell);
//FIXME          ufc.finite_elements[i]->evaluate_reference_basis_derivatives(j, 1, phi_grads[i][j][qp], q_points[qp], ref_cell);
          ufc.finite_elements[i]->evaluate_basis_derivatives(j, 1, phi_grads[i][j][qp], q_points[qp], ref_cell);
        }
      }
    }

    for (unsigned int i = 0; i<tensor_rank; ++i)
    {
//      std::cout << "i=" << i << std::endl;
      ref_phi_values[i].resize(n_test[i]);
      ref_phi_grads[i].resize(n_test[i]);
      for (unsigned int j = 0; j<n_test[i]; ++j)
      {
//        std::cout << "j=" << j << std::endl;
        ref_phi_values[i][j].resize(ref_points.size());
        ref_phi_grads[i][j].resize(ref_points.size());
        for (unsigned int qp = 0; qp<ref_points.size(); ++qp)
        {
          ref_phi_values[i][j][qp] = new real[1];
          ref_phi_grads[i][j][qp] = new real[dim];
          ufc.finite_elements[i]->evaluate_basis(j, ref_phi_values[i][j][qp], ref_points[qp], ref_cell);
          //FIXME:ufc.finite_elements[i]->evaluate_reference_basis_derivatives(j, 1, ref_phi_grads[i][j][qp], ref_points[qp], ref_cell);
          ufc.finite_elements[i]->evaluate_basis_derivatives(j, 1, ref_phi_grads[i][j][qp], ref_points[qp], ref_cell);
        }
      }
    }

    std::vector<std::vector<std::vector<double> > > Jinv_phi_i(q.size()); 
    std::vector<std::vector<std::vector<double> > > Jinv_phi_j(q.size());
    for (unsigned int ip = 0; ip < q.size(); ++ip)
    {
      const double detxW = det*weights[ip];
      Jinv_phi_i[ip].resize(n_test[0]);

      // Loop primary indices.
      for (unsigned int j = 0; j < n_test[0]; ++j)
      {
        Jinv_phi_i[ip][j].resize(dim);
        for(unsigned int d=0; d<dim; ++d)
          for(unsigned int e=0; e<dim; ++e)
            Jinv_phi_i[ip][j][d] += Jinv[e][d] * phi_grads[0][j][ip][e];
      }

      if(tensor_rank == 2)
      {
        Jinv_phi_j[ip].resize(n_test[1]);
        for (unsigned int k = 0; k < n_test[1]; ++k)
        {
          Jinv_phi_j[ip][k].resize(dim);
          for(unsigned int d=0; d<dim; ++d)
            for(unsigned int e=0; e<dim; ++e)
              Jinv_phi_j[ip][k][d] += Jinv[e][d] * phi_grads[1][k][ip][e];
        }
      }
    }

    // Loop quadrature points (tensor/monomial terms (0,))
    for (unsigned int ip = 0; ip < q.size(); ++ip)
    {
      const double Gip0 = det*weights[ip];
      
      // Loop primary indices.
      for (unsigned int j = 0; j < n_test[0]; ++j)
      {
//        if(std::abs(Jinv_phi_i[ip][j][0])<1.e-15)
//          Jinv_phi_i[ip][j][0] = 0.;
//        if(std::abs(Jinv_phi_i[ip][j][1])<1.e-15)
//          Jinv_phi_i[ip][j][1] = 0.;

        if(tensor_rank == 2)
        {
          for (unsigned int k = 0; k < n_test[1]; ++k)
          {
            real scalar_product = 0.;
            for(unsigned int d=0; d<dim; ++d)
              scalar_product += Jinv_phi_i[ip][j][d]*Jinv_phi_j[ip][k][d];

//            if(std::abs(scalar_product*Gip0)>1.e-15)
              ufc.A[j*n_test[1] + k] += scalar_product*Gip0;
          }// end loop over 'k'
        }
        else if(tensor_rank == 1)
        {
//          if(std::abs(phi_values[0][j][ip][0]*Gip0)>1.e-15)
          ufc.A[j] += phi_values[0][j][ip][0]*q_coefficients[0][ip][0]*Gip0;
        }
      }// end loop over 'j'
      
    }// end loop over 'ip

   for(unsigned int i = 0; i<real_points.size(); ++i)
      delete real_points[i];

    for(unsigned int i = 0; i<q_coefficients.size(); ++i)
      for(unsigned int q = 0; q<q_coefficients[i].size(); ++q)
        delete q_coefficients[i][q];

    for (unsigned int i = 0; i<tensor_rank; ++i)
      for (unsigned int j = 0; j<n_test[i]; ++j)
        for (unsigned int qp = 0; qp<q.size(); ++qp)
        {
          delete phi_values[i][j][qp];
          delete phi_grads[i][j][qp];
        }
  }
}
#endif
