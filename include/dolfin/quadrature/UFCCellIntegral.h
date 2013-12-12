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

  class UFCCellIntegral  //: public ufc::cell_integral
  {
  public:
    
    /// Constructor
    UFCCellIntegral();

    /// Constructor
//    UFCCellIntegral(unsigned int n);

    /// Constructor
//    UFCCellIntegral(unsigned int n, unsigned int m);

    /// Destructor
    virtual ~UFCCellIntegral();

    ///Tabulate the tensor for the contribution from a local cell
//    virtual void tabulate_tensor(double* A,
//        const double * const * w,
//        const ufc::cell& c) const;
    
    ///Tabulate the tensor for the contribution from a local cell
    virtual void tabulate_tensor(double* A,
        const double * const * w,
        const ufc::cell& c) const;
    
    ///Tabulate the tensor for the contribution from a local cell
    virtual void tabulate_tensor(const UFC& ufc, Cell& cell,
        const Array<Function*>& coefficients,
        const QuadratureRule& q, unsigned int** dofs) const;

    /// Return number of quadrature points
//    int size() const;

    /// Return quadrature point
//    real point(unsigned int i) const;

    /// Return quadrature weight
//    real weight(unsigned int i) const;

    /// Return sum of weights (length, area, volume)
//    real measure() const;
    
    /// Display quadrature data
//    virtual void disp() const;

  protected:
    
//    uint n;        // Number of quadrature points
//    real* points;  // Quadrature points
//    real* weights; // Quadrature weights
//    real m;        // Sum of weights
  };

  inline UFCCellIntegral::UFCCellIntegral()
  {
    std::cout << "new integral" << std::endl;
  }

//  UFCCellIntegral::UFCCellIntegral(unsigned int n)
//  {}

//  UFCCellIntegral::UFCCellIntegral(unsigned int n, unsigned int m)
//  {}

  inline UFCCellIntegral::~UFCCellIntegral()
  {}

  inline void UFCCellIntegral::tabulate_tensor(double* A,
        const double * const * w,
        const ufc::cell& c) const
  {
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
      
    // Compute determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;
      
    // Compute inverse of Jacobian
    const double Jinv_00 =  J_11 / detJ;
    const double Jinv_01 = -J_01 / detJ;
    const double Jinv_10 = -J_10 / detJ;
    const double Jinv_11 =  J_00 / detJ;
    
    // Set scale factor
    const double det = std::abs(detJ);
    
    
    // Reset values of the element tensor block
    for (unsigned int j = 0; j < 3; j++)
    {
      for (unsigned int k = 0; k < 3; k++)
      {
        A[j*3 + k] = 0;
      }// end loop over 'k'
    }// end loop over 'j'
    
    // Array of quadrature weights (tensor/monomial term 0)
    const static double W0[4] = {0.159020690871988, 0.0909793091280113, 0.159020690871988, 0.0909793091280113};
    // Array of quadrature weights (tensor/monomial term 1)
    const static double W1 = 0.5;
    
    const static double P_t1_p0_s0[1][2] = \
    {{-1, 1}};
    // Array of non-zero columns
    static const unsigned int nzc1[2] = {0, 2};
    // Array of non-zero columns
    static const unsigned int nzc0[2] = {0, 1};
    
    const static double P_t0_p1[4][3] = \
    {{0.666390246014701, 0.178558728263616, 0.155051025721682},
    {0.280019915499074, 0.0750311102226082, 0.644948974278318},
    {0.178558728263616, 0.666390246014701, 0.155051025721682},
    {0.0750311102226081, 0.280019915499074, 0.644948974278318}};
    
    // Number of operations to compute geometry constants = 25
    const double G0 = Jinv_00*Jinv_10*W1*det*w[1][0];
    const double G1 = Jinv_01*Jinv_11*W1*det*w[1][0];
    const double G2 = Jinv_00*Jinv_00*W1*det*w[1][0];
    const double G3 = Jinv_01*Jinv_01*W1*det*w[1][0];
    const double G4 = Jinv_10*Jinv_10*W1*det*w[1][0];
    const double G5 = Jinv_11*Jinv_11*W1*det*w[1][0];
    const double G6 = 1.0/w[0][0]*det;
    
    // Loop quadrature points (tensor/monomial terms (1,))
    // Number of operations to compute element tensor for following IP loop = 51
    // Only 1 integration point, omitting IP loop.
    
    // Number of operations to compute declarations = 3
    const double Gip0 = G0 + G1;
    const double Gip1 = G2 + G3;
    const double Gip2 = G4 + G5;
    
    // Loop primary indices.
    // Number of operations for primary indices = 48
    for (unsigned int j = 0; j < 2; j++)
    {
      for (unsigned int k = 0; k < 2; k++)
      {
        // Number of operations to compute entry = 3
        A[nzc0[j]*3 + nzc1[k]] += P_t1_p0_s0[0][j]*P_t1_p0_s0[0][k]*Gip0;
        // Number of operations to compute entry = 3
        A[nzc0[j]*3 + nzc0[k]] += P_t1_p0_s0[0][j]*P_t1_p0_s0[0][k]*Gip1;
        // Number of operations to compute entry = 3
        A[nzc1[j]*3 + nzc1[k]] += P_t1_p0_s0[0][j]*P_t1_p0_s0[0][k]*Gip2;
        // Number of operations to compute entry = 3
        A[nzc1[j]*3 + nzc0[k]] += P_t1_p0_s0[0][j]*P_t1_p0_s0[0][k]*Gip0;
      }// end loop over 'k'
    }// end loop over 'j'
    
    // Loop quadrature points (tensor/monomial terms (0,))
    // Number of operations to compute element tensor for following IP loop = 112
    for (unsigned int ip = 0; ip < 4; ip++)
    {
      
      // Number of operations to compute declarations = 1
      const double Gip0 = G6*W0[ip];
      
      // Loop primary indices.
      // Number of operations for primary indices = 27
      for (unsigned int j = 0; j < 3; j++)
      {
        for (unsigned int k = 0; k < 3; k++)
        {
          // Number of operations to compute entry = 3
          A[j*3 + k] += P_t0_p1[ip][j]*P_t0_p1[ip][k]*Gip0;
        }// end loop over 'k'
      }// end loop over 'j'
      
    }// end loop over 'i
  }


  inline void UFCCellIntegral::tabulate_tensor(const UFC& ufc, Cell& cell, 
      const Array<Function*>& coefficients,
      const QuadratureRule& q,
      unsigned int** dofs) const
  {
    // Extract vertex coordinates
    const double * const * x = ufc.cell.coordinates;
    
    UFCReferenceCell ref_cell(cell);
    const double * const * x_ref = ref_cell.coordinates;
    // Compute Jacobian of affine map from reference cell
//    const double J_00 = x[1][0] - x[0][0];
//    const double J_01 = x[2][0] - x[0][0];
//    const double J_10 = x[1][1] - x[0][1];
//    const double J_11 = x[2][1] - x[0][1];


    const unsigned int dim = ufc.mesh.topological_dimension;
    const unsigned int n_vertices = (dim == 1 ? 2 : (dim == 2 ? 3 : 4));
//    std::vector<std::vector<double> > x(n_vertices, std::vector<double>(dim,0.));

//    x[0][0] = 0;
//    x[0][1] = 0;
//    x[1][0] = 1;
//    x[1][1] = 0;
//    x[2][0] = 0;
//    x[2][1] = 1;

//    std::cout << "x0 = " << x[0][0] << "  " << x[0][1] <<  ",    x1 = " << x[1][0] << "  " << x[1][1]  << ",    x2 = " << x[2][0] << "  " << x[2][1] << std::endl;

//    std::cout << "x ref 0 = " << x_ref[0][0] << "  " << x_ref[0][1] << std::endl;
//    std::cout << "x ref 1 = " << x_ref[1][0] << "  " << x_ref[1][1] << std::endl;
//    std::cout << "x ref 2 = " << x_ref[2][0] << "  " << x_ref[2][1] << std::endl;

    std::vector<std::vector<double> > J (2, std::vector<double>(2,0.));
    J[0][0] = x[1][0] - x[0][0];
    J[0][1] = x[2][0] - x[0][0];
    J[1][0] = x[1][1] - x[0][1];
    J[1][1] = x[2][1] - x[0][1];

//    std::cout << "J_00 = " << J[0][0] << ",  J_01 = " << J[0][1] << ",  J_10 = " << J[1][0] << ",  J_11 = " << J[1][1] << std::endl;
    // Compute determinant of Jacobian
    double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
      
    // Compute inverse of Jacobian
    std::vector<std::vector<double> > Jinv (2, std::vector<double>(2,0.));
    Jinv[0][0] =  J[1][1] / detJ;
    Jinv[0][1] = -J[0][1] / detJ;
    Jinv[1][0] = -J[1][0] / detJ;
    Jinv[1][1] =  J[0][0] / detJ;
//    std::cout << "J_00 = " << Jinv[0][0] << std::endl;
//    std::cout << "J_01 = " << Jinv[0][1] << std::endl;
//    std::cout << "J_10 = " << Jinv[1][0] << std::endl;
//    std::cout << "J_11 = " << Jinv[1][1] << std::endl;
    
    // Set scale factor
    const double det = std::abs(detJ);
//    std::cout << "det=" << det << std::endl;
    
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
//      std::cout << "i= " << i << "  ";
        real_points[i] = new real[dim];
        for(unsigned int d=0; d<dim; ++d)
        {
          real_points[i][d] = 0.;
          for(unsigned int e=0; e<dim; ++e)
            real_points[i][d] += J[d][e] * q_points[i][e];
          real_points[i][d] += x[0][d];
        }
//        std::cout << real_points[i][0] << "  " << real_points[i][1] << std::endl;
    }

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
    
    std::vector<std::vector<real*> > q_coefficients (coefficients.size());
    std::vector<std::vector<real*> > ref_q_coefficients (coefficients.size());
    for(unsigned int i = 0; i<q_coefficients.size(); ++i)
    {
      q_coefficients[i].resize(q.size());
      for(unsigned int q = 0; q<q_coefficients[i].size(); ++q)
      {
        q_coefficients[i][q] = new real[1];
        coefficients[i]->eval(q_coefficients[i][q], real_points[q]);//, ref_cell);
//        coefficients[i]->eval(q_coefficients[i][q], q_points[q]);//, ref_cell);
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
//    std::cout << "coefs" << std::endl;
//    for(unsigned int i = 0; i<q_coefficients.size(); ++i)
//      for(unsigned int q = 0; q<q_coefficients[i].size(); ++q)
//        std::cout << q_coefficients[i][q][0] << "  ";
//    std::cout << std::endl;

//    std::cout << "ref points" << std::endl;
//    for (unsigned int qp = 0; qp<ref_points.size(); ++qp)
//      for(unsigned int d=0; d<dim; ++d)
//        std::cout << ref_points[qp][d] << "  ";
//    std::cout << std::endl;


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
//          std::cout << "phi = " << phi_values[i][j][qp][0] << std::endl;
//FIXME          ufc.finite_elements[i]->evaluate_reference_basis_derivatives(j, 1, phi_grads[i][j][qp], q_points[qp], ref_cell);
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
          //FIXME:ufc.finite_elements[i]->evaluate_reference_basis_derivatives(j, 1, ref_phi_grads[i][j][qp], ref_points[qp], ref_cell);
          ufc.finite_elements[i]->evaluate_basis_derivatives(j, 1, ref_phi_grads[i][j][qp], ref_points[qp], ref_cell);
//          std::cout << "grad phi(" << ref_points[qp][0] << "," << ref_points[qp][1] << ") = " 
//            << ref_phi_grads[i][j][qp][0] << "  " << ref_phi_grads[i][j][qp][1] << std::endl;
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
//      std::cout << "ip = " << ip << " w(ip) = " << weights[ip] << " det = " << det << std::endl;
      const double Gip0 = det*weights[ip];
      
      // Loop primary indices.
      for (unsigned int j = 0; j < n_test[0]; ++j)
      {
//        if(std::abs(Jinv_phi_i[ip][j][0])<1.e-15)
//          Jinv_phi_i[ip][j][0] = 0.;
//        if(std::abs(Jinv_phi_i[ip][j][1])<1.e-15)
//          Jinv_phi_i[ip][j][1] = 0.;

//        std::cout << "j = " << j << "  " //<< phi_grads[0][j][ip][0] << "  " << phi_grads[0][j][ip][1] << "  "
//          << Jinv_phi_i[ip][j][0] << "  " << Jinv_phi_i[ip][j][1] << std::endl;
        if(tensor_rank == 2)
        {
          for (unsigned int k = 0; k < n_test[1]; ++k)
          {
//            std::cout << "k = " << k << "  " << phi_grads[0][k][ip][0] << "  " << phi_grads[0][k][ip][1] << "  " 
//              << Jinv_phi_i[ip][k][0] << "  " << Jinv_phi_i[ip][k][1]<< std::endl;
            real scalar_product = 0.;
            for(unsigned int d=0; d<dim; ++d)
              scalar_product += Jinv_phi_i[ip][j][d]*Jinv_phi_j[ip][k][d];

//            if(std::abs(scalar_product*Gip0)>1.e-15)
              ufc.A[j*n_test[1] + k] += scalar_product*Gip0;
//              ufc.A[j*n_test[1] + k] += Jinv_phi_j[ip][k][0]*phi_values[0][j][ip][0]*Gip0;
//              ufc.A[j*n_test[1] + k] += Jinv_phi_i[ip][j][0]*phi_values[1][k][ip][0]*Gip0;
//              ufc.A[j*n_test[1] + k] += phi_values[0][j][ip][0]*phi_values[1][k][ip][0]*Gip0; //nicht exakt integrierbar
//            std::cout << "j = " << j << " k = " << k << "  " << scalar_product << std::endl;

          }// end loop over 'k'
        }
        else if(tensor_rank == 1)
        {
          ufc.A[j] += phi_values[0][j][ip][0]*q_coefficients[0][ip][0]*Gip0;
//          std::cout << "A(" <<  j << ")= " << ufc.A[j] << std::endl;
        }
//          if(std::abs(phi_values[0][j][ip][0]*Gip0)>1.e-15)
//            ufc.A[j] += Jinv_phi_i[ip][j][0]*ufc.w[0][j]*Gip0;
      }// end loop over 'j'
      
    }// end loop over 'ip

    for (unsigned int j = 0; j < n_test[0]; ++j)
    {
      if(tensor_rank == 2)
      {
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
      delete real_points[i];

    for(unsigned int i = 0; i<q_coefficients.size(); ++i)
      for(unsigned int q = 0; q<q_coefficients[i].size(); ++q)
        delete q_coefficients[i][q];

    for(unsigned int i = 0; i<ref_q_coefficients.size(); ++i)
      for(unsigned int q = 0; q<ref_q_coefficients[i].size(); ++q)
        delete ref_q_coefficients[i][q];

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
