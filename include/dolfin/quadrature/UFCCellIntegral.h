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
<<<<<<< HEAD
=======
//    std::vector<std::vector<double> > x(n_vertices, std::vector<double>(dim,0.));

//    x[0][0] = 0;
//    x[0][1] = 0;
//    x[1][0] = 1;
//    x[1][1] = 0;
//    x[2][0] = 0;
//    x[2][1] = 1;

<<<<<<< HEAD
//    std::cout << "x0 = " << x[0][0] << "  " << x[0][1] <<  ",    x1 = " << x[1][0] << "  " << x[1][1]  << ",    x2 = " << x[2][0] << "  " << x[2][1] << std::endl;

=======
//    std::cout << "x0 = " << x[0][0] << "  " << x[0][1] << std::endl;
//    std::cout << "x1 = " << x[1][0] << "  " << x[1][1] << std::endl;
//    std::cout << "x2 = " << x[2][0] << "  " << x[2][1] << std::endl;
//
>>>>>>> medium stage in adding quadrature rules
//    std::cout << "x ref 0 = " << x_ref[0][0] << "  " << x_ref[0][1] << std::endl;
//    std::cout << "x ref 1 = " << x_ref[1][0] << "  " << x_ref[1][1] << std::endl;
//    std::cout << "x ref 2 = " << x_ref[2][0] << "  " << x_ref[2][1] << std::endl;

>>>>>>> commit quadrature
    std::vector<std::vector<double> > J (2, std::vector<double>(2,0.));
    J[0][0] = x[1][0] - x[0][0];
    J[0][1] = x[2][0] - x[0][0];
    J[1][0] = x[1][1] - x[0][1];
    J[1][1] = x[2][1] - x[0][1];
<<<<<<< HEAD

=======
      
<<<<<<< HEAD
>>>>>>> commit quadrature
//    std::cout << "J_00 = " << J[0][0] << ",  J_01 = " << J[0][1] << ",  J_10 = " << J[1][0] << ",  J_11 = " << J[1][1] << std::endl;
=======
//    std::cout << "J_00 = " << J_00 << std::endl;
//    std::cout << "J_01 = " << J_01 << std::endl;
//    std::cout << "J_10 = " << J_10 << std::endl;
//    std::cout << "J_11 = " << J_11 << std::endl;
>>>>>>> medium stage in adding quadrature rules
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
<<<<<<< HEAD
=======
<<<<<<< HEAD
//    std::cout << "det=" << det << std::endl;
=======
    std::cout << "det=" << det << std::endl;
>>>>>>> medium stage in adding quadrature rules
>>>>>>> commit quadrature
    
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

<<<<<<< HEAD
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

<<<<<<< HEAD
=======
    /*
    std::cout << "q_points" << std::endl;
    for(unsigned int i = 0; i<q_points.size(); ++i)
    {
      for(unsigned int d=0; d<dim; ++d)
        std::cout << q_points[i][d] << "   ";
      std::cout << std::endl;
    }

    std::cout << "real points" << std::endl;
    for(unsigned int i = 0; i<real_points.size(); ++i)
    {
      for(unsigned int d=0; d<dim; ++d)
        std::cout << real_points[i][d] << "   ";
      std::cout << std::endl;
    }
    */
    
=======
>>>>>>> medium stage in adding quadrature rules
>>>>>>> commit quadrature
    std::vector<std::vector<real*> > q_coefficients (coefficients.size());
    for(unsigned int i = 0; i<q_coefficients.size(); ++i)
    {
      q_coefficients[i].resize(q.size());
      for(unsigned int q = 0; q<q_coefficients[i].size(); ++q)
      {
        q_coefficients[i][q] = new real[1];
<<<<<<< HEAD
        coefficients[i]->eval(q_coefficients[i][q], real_points[q]);
      }
    }

=======
<<<<<<< HEAD
        coefficients[i]->eval(q_coefficients[i][q], real_points[q]);//, ref_cell);
//        coefficients[i]->eval(q_coefficients[i][q], q_points[q]);//, ref_cell);
=======
        coefficients[i]->eval(q_coefficients[i][q], q_points[q]);//, ref_cell);
>>>>>>> medium stage in adding quadrature rules
      }
    }

    for(unsigned int i = 0; i<ref_q_coefficients.size(); ++i)
    {
      ref_q_coefficients[i].resize(ref_points.size());
      for(unsigned int q = 0; q<ref_q_coefficients[i].size(); ++q)
      {
        ref_q_coefficients[i][q] = new real[1];
        coefficients[i]->eval(ref_q_coefficients[i][q], ref_points[q]);//, ref_cell);
      }
    }
<<<<<<< HEAD
//    std::cout << "coefs" << std::endl;
//    for(unsigned int i = 0; i<q_coefficients.size(); ++i)
//      for(unsigned int q = 0; q<q_coefficients[i].size(); ++q)
//        std::cout << q_coefficients[i][q][0] << "  ";
//    std::cout << std::endl;
=======
    std::cout << "coefs" << std::endl;
    for(unsigned int i = 0; i<q_coefficients.size(); ++i)
      for(unsigned int q = 0; q<q_coefficients[i].size(); ++q)
        std::cout << q_coefficients[i][q][0] << "  ";
    std::cout << std::endl;
>>>>>>> medium stage in adding quadrature rules

//    std::cout << "ref points" << std::endl;
//    for (unsigned int qp = 0; qp<ref_points.size(); ++qp)
//      for(unsigned int d=0; d<dim; ++d)
//        std::cout << ref_points[qp][d] << "  ";
//    std::cout << std::endl;


>>>>>>> commit quadrature
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
<<<<<<< HEAD
          ufc.finite_elements[i]->evaluate_reference_basis_derivatives(j, 1, phi_grads[i][j][qp], q_points[qp], ref_cell);
=======
//          std::cout << "phi = " << phi_values[i][j][qp][0] << std::endl;
<<<<<<< HEAD
//FIXME          ufc.finite_elements[i]->evaluate_reference_basis_derivatives(j, 1, phi_grads[i][j][qp], q_points[qp], ref_cell);
=======
>>>>>>> medium stage in adding quadrature rules
          ufc.finite_elements[i]->evaluate_basis_derivatives(j, 1, phi_grads[i][j][qp], q_points[qp], ref_cell);
//          std::cout << "grad phi = " << phi_grads[i][j][qp][0] << "  " << phi_grads[i][j][qp][1] << std::endl;
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
//          std::cout << "phi(" << ref_points[qp][0] << "," << ref_points[qp][1] << ") = " << ref_phi_values[i][j][qp][0] << std::endl;
<<<<<<< HEAD
          //FIXME:ufc.finite_elements[i]->evaluate_reference_basis_derivatives(j, 1, ref_phi_grads[i][j][qp], ref_points[qp], ref_cell);
=======
>>>>>>> medium stage in adding quadrature rules
          ufc.finite_elements[i]->evaluate_basis_derivatives(j, 1, ref_phi_grads[i][j][qp], ref_points[qp], ref_cell);
//          std::cout << "grad phi(" << ref_points[qp][0] << "," << ref_points[qp][1] << ") = " 
//            << ref_phi_grads[i][j][qp][0] << "  " << ref_phi_grads[i][j][qp][1] << std::endl;
>>>>>>> commit quadrature
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

<<<<<<< HEAD
   for(unsigned int i = 0; i<real_points.size(); ++i)
=======
    for (unsigned int j = 0; j < n_test[0]; ++j)
    {
      if(tensor_rank == 2)
      {
<<<<<<< HEAD
//        for (unsigned int k = 0; k < n_test[1]; ++k)
//          std::cout << "A(" <<  j*n_test[1] + k << ")= " << ufc.A[j*n_test[1] + k] << std::endl;
      }
      else if(tensor_rank == 1)
      {
//        std::cout << "A(" <<  j << ")= " << ufc.A[j] << std::endl;
//        std::cout << "w[" <<  j << "]= " << ufc.w[0][j] << std::endl;
      }
    } 

    for(unsigned int i = 0; i<real_points.size(); ++i)
>>>>>>> commit quadrature
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
=======
        for (unsigned int k = 0; k < n_test[1]; ++k)
          std::cout << "A(" <<  j*n_test[1] + k << ")= " << ufc.A[j*n_test[1] + k] << std::endl;
      }
      else if(tensor_rank == 1)
      {
        std::cout << "A(" <<  j << ")= " << ufc.A[j] << std::endl;
//        std::cout << "w[" <<  j << "]= " << ufc.w[0][j] << std::endl;
      }
    } 
  }

>>>>>>> medium stage in adding quadrature rules
}
#endif
