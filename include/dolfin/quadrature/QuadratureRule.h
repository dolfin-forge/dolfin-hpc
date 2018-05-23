// Copyright (C) 2003-2006 Bärbel Janssen.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __DOLFIN_QUADRATURE_RULE_H
#define __DOLFIN_QUADRATURE_RULE_H

#include <dolfin/common/types.h>

namespace dolfin
{

class QuadratureRule
{
public:

  /// Return number of quadrature points
  unsigned int size() const;

  /// Return quadrature point
  const real* point(unsigned int i) const;

  /// Return quadrature weight
  real weight(unsigned int i) const;

  /// Return quadrature points
  const std::vector<real*> & get_points() const;

  /// Return points of reference element
  const std::vector<real*> & get_reference_points() const;

  /// Return quadrature weight
  const std::vector<real> & get_weights() const;

  /// Return sum of weights (length, area, volume)
  real measure() const;

  /// Display quadrature data
  void disp() const;

protected:

  /// Constructor
  //QuadratureRule(unsigned int dim, unsigned int n);

  /// Destructor
  virtual ~QuadratureRule();

//      uint n;        // Number of quadrature points
  uint dim;        // Dimension of quadrature rule
  std::vector<real*> points;  // Quadrature points
  std::vector<real*> reference_points;  // Points of reference element
  std::vector<real> weights;  // Quadrature weights
  real m;        // Sum of weights

  void fill_reference_points();

};

//-----------------------------------------------------------------------------
inline QuadratureRule::~QuadratureRule()
{
  for (unsigned int i = 0; i < points.size(); ++i)
    delete points[i];
  for (unsigned int i = 0; i < reference_points.size(); ++i)
    delete reference_points[i];
}

//-----------------------------------------------------------------------------
inline unsigned int QuadratureRule::size() const
{
  return weights.size();
}

//-----------------------------------------------------------------------------
inline const real* QuadratureRule::point(unsigned int i) const
{
  if (i > points.size()) error("Invalid quadrature point");
  return points[i];
}

//-----------------------------------------------------------------------------
inline real QuadratureRule::weight(unsigned int i) const
{
  if (i > weights.size()) error("Invalid quadrature point");
  return weights[i];
}

//-----------------------------------------------------------------------------
inline const std::vector<real*> & QuadratureRule::get_points() const
{
  return points;
}

//-----------------------------------------------------------------------------
inline const std::vector<real*> & QuadratureRule::get_reference_points() const
{
  return reference_points;
}

//-----------------------------------------------------------------------------
inline const std::vector<real> & QuadratureRule::get_weights() const
{
  return weights;
}

//-----------------------------------------------------------------------------
inline real QuadratureRule::measure() const
{
  return m;
}

//-----------------------------------------------------------------------------
inline void QuadratureRule::disp() const
{
  cout << "QuadratureRule\n";
  cout << "--------------\n";
  begin("");
  cout << "Size             : " << uint(weights.size()) << "\n";
  cout << "Scaling factor   : " << m << "\n";
  cout << "Weights          : \n" << "\n";
  for (unsigned int i = 0; i < weights.size(); ++i)
  {
    cout << "\t" << weights[i];
  }
  cout << "\n";
  cout << "\n";
  cout << "Points           : \n" << "\n";
  for (unsigned int i = 0; i < points.size(); ++i)
  {
    for (unsigned int d = 0; d < dim; ++d)
    {
      cout << "\t" << points[i][d];
    }
    cout << "\n";
  }
  cout << "\n";
  cout << "Reference points : \n" << "\n";
  for (unsigned int i = 0; i < points.size(); ++i)
  {
    for (unsigned int d = 0; d < dim; ++d)
    {
      cout << "\t" << reference_points[i][d];
    }
    cout << "\n";
  }
  cout << "\n";
  end();
}

inline void QuadratureRule::fill_reference_points()
{
  reference_points.resize(dim + 1);
  switch (dim)
    {
    case 1:
      {
        for (unsigned int i = 0; i < reference_points.size(); ++i)
          reference_points[i] = new real[dim];

        for (unsigned int i = 0; i < reference_points.size(); ++i)
          reference_points[i][0] = i;

        break;
      }
    case 2:
      {
        for (unsigned int i = 0; i < reference_points.size(); ++i)
          reference_points[i] = new real[dim];

        unsigned int k = 0;
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            if (i + 1 + j >= reference_points.size()) break;
            else
            {
              reference_points[k][0] = i;
              reference_points[k][1] = j;
              k++;
            }
        break;
      }
    case 3:
      {
        for (unsigned int i = 0; i < reference_points.size(); ++i)
          reference_points[i] = new real[dim];

        unsigned int k = 0;
        for (unsigned int i = 0; i < dim; ++i)
          for (unsigned int j = 0; j < dim; ++j)
            for (unsigned int h = 0; h < dim; ++h)
              if (i + j + h + 1 >= reference_points.size()) break;
              else
              {
                reference_points[k][0] = i;
                reference_points[k][1] = j;
                reference_points[k][2] = h;
                k++;
              }
        break;
      }
    default:
      error("Reference points are available only from dimension 1 to 3.");
      break;
    }
}
}

#endif
