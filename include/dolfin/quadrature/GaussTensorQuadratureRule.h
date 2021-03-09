// Copyright (C) 2003-2005 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/config/dolfin_config.h>


#ifndef __DOLFIN_GAUSS_TENSOR_QUADRATURE_RULE_H
#define __DOLFIN_GAUSS_TENSOR_QUADRATURE_RULE_H

#include "QuadratureRule.h"

#include <vector>

namespace dolfin
{

  /// Gaussian-type quadrature rule on the reference cell.
  ///
  /// Points and weights are computed using one dimensional
  /// GaussQuadrature.

  class GaussTensorQuadratureRule : public QuadratureRule
  {
    public:

      GaussTensorQuadratureRule(unsigned int dim, unsigned int n);

      GaussTensorQuadratureRule(unsigned int dim, std::vector<unsigned int>& n);
  };

  inline GaussTensorQuadratureRule::GaussTensorQuadratureRule(unsigned int dim, unsigned int degree)
  {
    this->dim = dim;
    fill_reference_points();
    m = 0.;
    //create one dimensional quadratures
    std::vector<GaussQuadrature*> quadratures(1);
    quadratures[0] = new GaussQuadrature(degree);

    unsigned int n = 1;
    for(unsigned int d=0; d<dim; ++d)
      n *= degree;

    points.resize(n);
    for(unsigned int i=0; i<n; ++i)
      points[i] = new real[dim];
    weights.resize(n);

    unsigned int size_1d = quadratures[0]->size();
    switch(dim)
    {
      case 1:
        {
          for(unsigned int i=0; i<size_1d; ++i)
          {
            points[i][0] = quadratures[0]->point(i);
            weights[i] = quadratures[0]->weight(i);
          }
          break;
        }
      case 2:
        {
          for(unsigned int i=0; i<size_1d; ++i)
            for(unsigned int j=0; j<size_1d; ++j)
            {
              points[i*size_1d+j][0] = .5*(1+quadratures[0]->point(i));
              points[i*size_1d+j][1] = .25*(1-quadratures[0]->point(i))
                *(1+quadratures[0]->point(j));
              weights[i*size_1d+j] = 0.125*quadratures[0]->weight(i)
                *quadratures[0]->weight(j)*(1-quadratures[0]->point(i));
            }
          break;
        }
      case 3:
        {
          for(unsigned int i=0; i<size_1d; ++i)
            for(unsigned int j=0; j<size_1d; ++j)
              for(unsigned int k=0; k<size_1d; ++k)
              {
                const unsigned int index = (i*size_1d+j)*size_1d+k;
                points[index][0] = .5*(1+quadratures[0]->point(i));
                points[index][1] = .25*(1-quadratures[0]->point(i))
                  *(1+quadratures[0]->point(j));
                points[index][2] = .125*(1-quadratures[0]->point(i))
                  *(1-quadratures[0]->point(j))
                  *(1+quadratures[0]->point(k));
                weights[index] = 0.015625*quadratures[0]->weight(i)
                  *quadratures[0]->weight(j)*quadratures[0]->weight(k)
                  *(1-quadratures[0]->point(i))*(1-quadratures[0]->point(i))
                  *(1-quadratures[0]->point(j));
              }
          std::cout << std::endl;
          break;
        }
      default:
        error("This case is not covered.");
    }

    for(unsigned int i=0; i<weights.size(); ++i)
      m += weights[i];

    delete quadratures[0];
  }

  inline GaussTensorQuadratureRule::GaussTensorQuadratureRule(unsigned int dim, std::vector<unsigned int>& degree)
  {
    if(degree.size() != dim)
      error("The size of the vector does not match the dimension");

    this->dim = dim;
    fill_reference_points();

    //create one dimensional quadratures
    std::vector<GaussQuadrature*> quadratures(dim);
    unsigned int n = 1;
    for(unsigned int i = 0; i<dim; ++i)
    {
      n *= degree[i];
      quadratures[i] = new GaussQuadrature(degree[i]);
    }

    points.resize(n);
    for(unsigned int i=0; i<n; ++i)
      points[i] = new real[dim];
    weights.resize(n);

    switch(dim)
    {
      case 1:
        {
          for(unsigned int i=0; i<quadratures[0]->size(); ++i)
          {
            points[i][0] = quadratures[0]->point(i);
            weights[i] = quadratures[0]->weight(i);
          }
          break;
        }
      case 2:
        {
          for(unsigned int i=0; i<quadratures[0]->size(); ++i)
            for(unsigned int j=0; j<quadratures[1]->size(); ++j)
            {
              points[i*quadratures[0]->size()+j][0] = .5*(1+quadratures[0]->point(i));
              points[i*quadratures[0]->size()+j][1] = .25*(1-quadratures[0]->point(i))
                *(1+quadratures[0]->point(j));
              weights[i*quadratures[0]->size()+j] = 0.125*quadratures[0]->weight(i)
                *quadratures[1]->weight(j)*(1-quadratures[0]->point(i));
            }
          break;
        }
      case 3:
        {
          for(unsigned int i=0; i<quadratures[0]->size(); ++i)
            for(unsigned int j=0; j<quadratures[1]->size(); ++j)
              for(unsigned int k=0; k<quadratures[2]->size(); ++k)
              {
                const unsigned int index = i*quadratures[0]->size() +
                  j*quadratures[1]->size() + k;
                points[index][0] = .5*(1+quadratures[0]->point(i));
                points[index][1] = .25*(1-quadratures[1]->point(i))
                  *(1+quadratures[0]->point(j));
                points[index][2] = .125*(1-quadratures[0]->point(i))
                  *(1-quadratures[1]->point(j))
                  *(1+quadratures[2]->point(k));
                weights[index] = .015625*quadratures[0]->weight(i)
                  *quadratures[1]->weight(j)*quadratures[2]->weight(k)
                  *(1-quadratures[0]->point(i))*(1-quadratures[0]->point(i))
                  *(1-quadratures[1]->point(j));
              }
          break;
        }
      default:
        error("This case is not covered.");
    }

    for(unsigned int i=0; i<weights.size(); ++i)
      m += weights[i];

    for(unsigned int i=0; i<quadratures.size(); ++i)
      delete quadratures[i];
  }

}
#endif
