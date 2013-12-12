// Copyright (C) 2003-2005 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#include <dolfin/config/dolfin_config.h>


#ifndef __GAUSS_QUADRATURE_RULE_H
#define __GAUSS_QUADRATURE_RULE_H

//#include "Quadrature.h"

namespace dolfin 
{

  /// Gaussian-type quadrature rule on the reference cell.
  /// Points and weights are hardcoded from literature.

  class GaussQuadratureRule : public QuadratureRule
  {
    public:

      GaussQuadratureRule(unsigned int dim, unsigned int n);
  };

  inline GaussQuadratureRule::GaussQuadratureRule(unsigned int dim, unsigned int degree) 
  {
    this->dim = dim;
    fill_reference_points();
    m = 0.;
    switch (dim)
    {
      case 1:
        {
          GaussQuadrature* q = new GaussQuadrature(degree);
          unsigned int n = q->size();
          points.resize(n);
          weights.resize(n);
          for(unsigned int i=0; i<n; ++i)
          {
            points[i] = new real[1];
            points[i][0] = q->point(i); 
            weights[i] = q->weight(i); 
          }
          delete q;
          break;
        }
      case 2:
        {
          switch (degree)
          {
            case 1:
              {
                unsigned int n = 1;
                points.resize(n);
                weights.resize(n);
                points[0] = new real[dim];

                points[0][0] =  0.333333333333;
                points[0][1] =  0.333333333333;
                weights[0] =    0.500000000000;      
                break;
              }
            case 2:
              {
                unsigned int n = 3;
                points.resize(n);
                weights.resize(n);
                for(unsigned int i=0; i<n; ++i)
                  points[i] = new real[dim];

                points[0][0] =  0.666666666666666667;
                points[0][1] =  0.166666666666666667;
                points[1][0] =  0.166666666666666667;
                points[1][1] =  0.166666666666666667;
                points[2][0] =  0.166666666666666667;
                points[2][1] =  0.666666666666666667;

                weights[0] =    0.166666666666666666;      
                weights[1] =    0.166666666666666666;      
                weights[2] =    0.166666666666666666;      
                break;
              }
            case 3:
              {
                unsigned int n = 4;
                points.resize(n);
                weights.resize(n);
                for(unsigned int i=0; i<n; ++i)
                  points[i] = new real[dim];

                points[0][0] = 0.3333333333333333333;
                points[0][1] = 0.3333333333333333333;
                points[1][0] = 0.2000000000000000000;
                points[1][1] = 0.2000000000000000000;
                points[2][0] = 0.6000000000000000000;
                points[2][1] = 0.2000000000000000000;
                points[3][0] = 0.2000000000000000000;
                points[3][1] = 0.6000000000000000000;

                weights[0] =-0.2812500000000000000;
                weights[1] = 0.2604166666666666666;
                weights[2] = 0.2604166666666666666;
                weights[3] = 0.2604166666666666666;

                /*
                unsigned int n = 6;
                points.resize(n);
                weights.resize(n);
                for(unsigned int i=0; i<n; ++i)
                  points[i] = new real[dim];


               points[0][0] = 0.659027622374092; 
               points[0][1] = 0.231933368553031;
               points[1][0] = 0.659027622374092; 
               points[1][1] = 0.109039009072877;
               points[2][0] = 0.231933368553031; 
               points[2][1] = 0.659027622374092;
               points[3][0] = 0.231933368553031; 
               points[3][1] = 0.109039009072877;
               points[4][0] = 0.109039009072877; 
               points[4][1] = 0.659027622374092;
               points[5][0] = 0.109039009072877; 
               points[5][1] = 0.231933368553031;

               weights[0] = 1./12.;
               weights[1] = 1./12.;
               weights[2] = 1./12.;
               weights[3] = 1./12.;
               weights[4] = 1./12.;
               weights[5] = 1./12.;
                */
               break;
              }
            case 4:
              {
                unsigned int n = 6;
                points.resize(n);
                weights.resize(n);
                for(unsigned int i=0; i<n; ++i)
                  points[i] = new real[dim];

                points[0][0] = 0.108103018168070; 
                points[0][1] = 0.445948490915965;
                points[1][0] = 0.445948490915965; 
                points[1][1] = 0.445948490915965;
                points[2][0] = 0.445948490915965; 
                points[2][1] = 0.108103018168070;
                points[3][0] = 0.816847572980458; 
                points[3][1] = 0.091576213509771;
                points[4][0] = 0.091576213509771; 
                points[4][1] = 0.091576213509771;
                points[5][0] = 0.091576213509771; 
                points[5][1] = 0.816847572980458;

                weights[0] = 0.111690794839005; 
                weights[1] = 0.111690794839005; 
                weights[2] = 0.111690794839005; 
                weights[3] = 0.054975871827661; 
                weights[4] = 0.054975871827661; 
                weights[5] = 0.054975871827661; 
                break;
              }
            case 5:
              {
                unsigned int n = 7;
                points.resize(n);
                weights.resize(n);
                for(unsigned int i=0; i<n; ++i)
                  points[i] = new real[dim];

                points[0][0] = 0.33333333333333333333; 
                points[0][1] = 0.33333333333333333333;
                points[1][0] = 0.10128650732345633880;
                points[1][1] = 0.10128650732345633880;
                points[2][0] = 0.79742698535308732240;
                points[2][1] = 0.10128650732345633880;
                points[3][0] = 0.10128650732345633880;
                points[3][1] = 0.79742698535308732240;
                points[4][0] = 0.47014206410511508977;
                points[4][1] = 0.47014206410511508977;
                points[5][0] = 0.05971587178976982046;
                points[5][1] = 0.47014206410511508977;
                points[6][0] = 0.47014206410511508977;
                points[6][1] = 0.05971587178976982046;

                weights[0] = 0.11250000000000000000; 
                weights[1] = 0.06296959027241357629;
                weights[2] = 0.06296959027241357629;
                weights[3] = 0.06296959027241357629;
                weights[4] = 0.06619707639425309036;
                weights[5] = 0.06619707639425309036;
                weights[6] = 0.06619707639425309036;
                break;
              }
            case 6:
              {
                unsigned int n = 12;
                points.resize(n);
                weights.resize(n);
                for(unsigned int i=0; i<n; ++i)
                  points[i] = new real[dim];

                points[0][0]  = 0.501426509658180; 
                points[0][1]  = 0.249286745170910;
                points[1][0]  = 0.249286745170910;
                points[1][1]  = 0.249286745170910;
                points[2][0]  = 0.249286745170910;
                points[2][1]  = 0.501426509658180;
                points[3][0]  = 0.873821971016996;
                points[3][1]  = 0.063089014491502;
                points[4][0]  = 0.063089014491502;
                points[4][1]  = 0.063089014491502;
                points[5][0]  = 0.063089014491502;
                points[5][1]  = 0.873821971016996;
                points[6][0]  = 0.053145049844816;
                points[6][1]  = 0.310352451033785;
                points[7][0]  = 0.310352451033785;
                points[7][1]  = 0.053145049844816;
                points[8][0]  = 0.053145049844816;
                points[8][1]  = 0.636502499121399;
                points[9][0]  = 0.636502499121399;
                points[9][1]  = 0.053145049844816;
                points[10][0] = 0.636502499121399;
                points[10][1] = 0.310352451033785;
                points[11][0] = 0.310352451033785;
                points[11][1] = 0.636502499121399;

                weights[0]  = 0.058393137863189; 
                weights[1]  = 0.058393137863189;
                weights[2]  = 0.058393137863189;
                weights[3]  = 0.025422453185103;
                weights[4]  = 0.025422453185103;
                weights[5]  = 0.025422453185103;
                weights[6]  = 0.041425537809187;
                weights[7]  = 0.041425537809187;
                weights[8]  = 0.041425537809187;
                weights[9]  = 0.041425537809187;
                weights[10] = 0.041425537809187;
                weights[11] = 0.041425537809187;
                break;
              }
            case 7:
              {
                unsigned int n = 13;
                points.resize(n);
                weights.resize(n);
                for(unsigned int i=0; i<n; ++i)
                  points[i] = new real[dim];

                points[0][0]  = 0.333333333333; 
                points[0][1]  = 0.333333333333;
                points[1][0]  = 0.479308067841;
                points[1][1]  = 0.260345966079;
                points[2][0]  = 0.260345966079;
                points[2][1]  = 0.260345966079;
                points[3][0]  = 0.260345966079;
                points[3][1]  = 0.479308067841;
                points[4][0]  = 0.869739794195;
                points[4][1]  = 0.065130102902;
                points[5][0]  = 0.065130102902;
                points[5][1]  = 0.065130102902;
                points[6][0]  = 0.065130102902;
                points[6][1]  = 0.869739794195;
                points[7][0]  = 0.048690315425;
                points[7][1]  = 0.312865496004;
                points[8][0]  = 0.312865496004;
                points[8][1]  = 0.048690315425;
                points[9][0]  = 0.048690315425;
                points[9][1]  = 0.638444188569;
                points[10][0] = 0.638444188569;
                points[10][1] = 0.048690315425;
                points[11][0] = 0.638444188569;
                points[11][1] = 0.312865496004;
                points[12][0] = 0.312865496004;
                points[12][1] = 0.638444188569;

                weights[0]  = -0.074785022233; 
                weights[1]  =  0.087807628716;
                weights[2]  =  0.087807628716;
                weights[3]  =  0.087807628716;
                weights[4]  =  0.026673617804;
                weights[5]  =  0.026673617804;
                weights[6]  =  0.026673617804;
                weights[7]  =  0.038556880445;
                weights[8]  =  0.038556880445;
                weights[9]  =  0.038556880445;
                weights[10] =  0.038556880445;
                weights[11] =  0.038556880445;
                weights[12] =  0.038556880445;
                break;
              }
            default:
              error("Gauss quadrature only implemented up to degree 7 in 2d.");
          }   
          break;
        }
      case 3:
        {
          switch (degree)
          {
            case 1:
              {
                unsigned int n = 1;
                points.resize(n);
                weights.resize(n);
                points[0] = new real[dim];

                points[0][0] = 0.250000000000;
                points[0][1] = 0.250000000000;
                points[0][2] = 0.250000000000;
                weights[0] =   0.166666666667;
                break;
              }
            case 2:
              {
                unsigned int n = 4;
                points.resize(n);
                weights.resize(n);
                for(unsigned int i=0; i<n; ++i)
                  points[i] = new real[dim];

                points[0][0] = 0.13819660112501051518;
                points[0][1] = 0.13819660112501051518;
                points[0][2] = 0.13819660112501051518;
                points[1][0] = 0.58541019662496845446;
                points[1][1] = 0.13819660112501051518;
                points[1][2] = 0.13819660112501051518;
                points[2][0] = 0.13819660112501051518;
                points[2][1] = 0.58541019662496845446;
                points[2][2] = 0.13819660112501051518;
                points[3][0] = 0.13819660112501051518;
                points[3][1] = 0.13819660112501051518;
                points[3][2] = 0.58541019662496845446;

                weights[0] = 0.04166666666666666666;
                weights[1] = 0.04166666666666666666;
                weights[2] = 0.04166666666666666666;
                weights[3] = 0.04166666666666666666;
                break;
              }
            case 3:
              {
                unsigned int n = 5;
                points.resize(n);
                weights.resize(n);
                for(unsigned int i=0; i<n; ++i)
                  points[i] = new real[dim];

                points[0][0] = 0.25;
                points[0][1] = 0.25;
                points[0][2] = 0.25;
                points[1][0] = 0.1666666666666667;
                points[1][1] = 0.1666666666666667;
                points[1][2] = 0.1666666666666667;
                points[2][0] = 0.5000000000000000;
                points[2][1] = 0.1666666666666667;
                points[2][2] = 0.1666666666666667;
                points[3][0] = 0.1666666666666667;
                points[3][1] = 0.5000000000000000;
                points[3][2] = 0.1666666666666667;
                points[4][0] = 0.1666666666666667;
                points[4][1] = 0.1666666666666667;
                points[4][2] = 0.5000000000000000;

                weights[0] = -0.13333333333333333333;
                weights[1] =  0.07500000000000000000;
                weights[2] =  0.07500000000000000000;
                weights[3] =  0.07500000000000000000;
                weights[4] =  0.07500000000000000000;
                break;
              }
            case 4:
              {
                unsigned int n = 11;
                points.resize(n);
                weights.resize(n);
                for(unsigned int i=0; i<n; ++i)
                  points[i] = new real[dim];

                points[0][0]  = 0.2500000000000000;
                points[0][1]  = 0.2500000000000000;
                points[0][2]  = 0.2500000000000000;
                points[1][0]  = 0.7857142857142857;
                points[1][1]  = 0.0714285714285714;
                points[1][2]  = 0.0714285714285714;
                points[2][0]  = 0.0714285714285714;
                points[2][1]  = 0.0714285714285714;
                points[2][2]  = 0.0714285714285714;
                points[3][0]  = 0.0714285714285714;
                points[3][1]  = 0.0714285714285714;
                points[3][2]  = 0.7857142857142857;
                points[4][0]  = 0.0714285714285714;
                points[4][1]  = 0.7857142857142857;
                points[4][2]  = 0.0714285714285714;
                points[5][0]  = 0.1005964238332008;
                points[5][1]  = 0.3994035761667992;
                points[5][2]  = 0.3994035761667992;
                points[6][0]  = 0.3994035761667992;
                points[6][1]  = 0.1005964238332008;
                points[6][2]  = 0.3994035761667992;
                points[7][0]  = 0.3994035761667992;
                points[7][1]  = 0.3994035761667992;
                points[7][2]  = 0.1005964238332008;
                points[8][0]  = 0.3994035761667992;
                points[8][1]  = 0.1005964238332008;
                points[8][2]  = 0.1005964238332008;
                points[9][0]  = 0.1005964238332008;
                points[9][1]  = 0.3994035761667992;
                points[9][2]  = 0.1005964238332008;
                points[10][0] = 0.1005964238332008;
                points[10][1] = 0.1005964238332008;
                points[10][2] = 0.3994035761667992;

                weights[0]  = -0.01315555555555555;
                weights[1]  =  0.00762222222222222;
                weights[2]  =  0.00762222222222222;
                weights[3]  =  0.00762222222222222;
                weights[4]  =  0.00762222222222222;
                weights[5]  =  0.02488888888888888;
                weights[6]  =  0.02488888888888888;
                weights[7]  =  0.02488888888888888;
                weights[8]  =  0.02488888888888888;
                weights[9]  =  0.02488888888888888;
                weights[10] =  0.02488888888888888;
                break;
              }
            case 5:
              {
                unsigned int n = 15;
                points.resize(n);
                weights.resize(n);
                for(unsigned int i=0; i<n; ++i)
                  points[i] = new real[dim];

                points[0][0]   = 0.25000000000000000000;
                points[0][1]   = 0.25000000000000000000;
                points[0][2]   = 0.25000000000000000000;
                points[1][0]   = 0.09197107805272303278; 
                points[1][1]   = 0.09197107805272303278;
                points[1][2]   = 0.09197107805272303278;
                points[2][0]   = 0.09197107805272303278;
                points[2][1]   = 0.09197107805272303278;
                points[2][2]   = 0.81605784389455393443;
                points[3][0]   = 0.09197107805272303278;
                points[3][1]   = 0.81605784389455393443;
                points[3][2]   = 0.09197107805272303278;
                points[4][0]   = 0.81605784389455393443;
                points[4][1]   = 0.09197107805272303278;
                points[4][2]   = 0.09197107805272303278;
                points[5][0]   = 0.31979362782962990838;
                points[5][1]   = 0.31979362782962990838;
                points[5][2]   = 0.31979362782962990838;
                points[6][0]   = 0.31979362782962990838;
                points[6][1]   = 0.31979362782962990838;
                points[6][2]   = 0.36041274434074018323;
                points[7][0]   = 0.31979362782962990838;
                points[7][1]   = 0.36041274434074018323;
                points[7][2]   = 0.31979362782962990838;
                points[8][0]   = 0.36041274434074018323;
                points[8][1]   = 0.31979362782962990838;
                points[8][2]   = 0.31979362782962990838;
                points[9][0]   = 0.13819660112501051518;
                points[9][1]   = 0.13819660112501051518;
                points[9][2]   = 0.36180339887498948482;
                points[10][0]  = 0.13819660112501051518;
                points[10][1]  = 0.36180339887498948482;
                points[10][2]  = 0.13819660112501051518;
                points[11][0]  = 0.36180339887498948482;
                points[11][1]  = 0.13819660112501051518;
                points[11][2]  = 0.13819660112501051518;
                points[12][0]  = 0.13819660112501051518;
                points[12][1]  = 0.36180339887498948482;
                points[12][2]  = 0.36180339887498948482;
                points[13][0]  = 0.36180339887498948482;
                points[13][1]  = 0.13819660112501051518;
                points[13][2]  = 0.36180339887498948482;
                points[14][0]  = 0.36180339887498948482;
                points[14][1]  = 0.36180339887498948482;
                points[14][2]  = 0.13819660112501051518;

                weights[0]  = 0.01975308641975308641;
                weights[1]  = 0.01198951396316977000; 
                weights[2]  = 0.01198951396316977000;
                weights[3]  = 0.01198951396316977000;
                weights[4]  = 0.01198951396316977000;
                weights[5]  = 0.01151136787104539754;
                weights[6]  = 0.01151136787104539754;
                weights[7]  = 0.01151136787104539754;
                weights[8]  = 0.01151136787104539754;
                weights[9]  = 0.00881834215167548500; 
                weights[10] = 0.00881834215167548500;
                weights[11] = 0.00881834215167548500;
                weights[12] = 0.00881834215167548500;
                weights[13] = 0.00881834215167548500;
                weights[14] = 0.00881834215167548500;
                break;
              }
            default:
              error("Gauss quadrature only implemented up to degree 5 in 3d.");
          }
          break;
        }
      default:
        error("Gauss quadrature only implemented in dimensions 1-3.");
    }

    for(unsigned int i=0; i<weights.size(); ++i)
      m += weights[i];
  }

}
#endif
