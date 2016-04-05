//
//
//

#include <dolfin/mesh/MeshQuality.h>

#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/Vertex.h>

#include <cmath>

namespace dolfin
{

//-----------------------------------------------------------------------------
MeshQuality::MeshQuality(Mesh& mesh) :
    MeshDependent(mesh),
    mu_min(0.0),
    mu_max(0.0),
    mu_avg(0.0),
    h_min(0.0),
    h_max(0.0),
    h_avg(0.0),
    vol_min(0.0),
    vol_max(0.0),
    orientation_(mesh, mesh.topology().dim()),
    mapping_(mesh),
    bbox_min_(),
    bbox_max_()
{
  // Initialize orientation
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    orientation_.set(c->index(), (uint) c->orientation());
  }

  //
  this->compute();
}

//-----------------------------------------------------------------------------
bool MeshQuality::is_inverted(uint& first)
{
  for (CellIterator c(mesh()); !c.end(); ++c)
  {
    if (orientation_.get(c->index()) != c->orientation())
    {
      first = c->index();
      return true;
    }
  }
  return false;
}

//-----------------------------------------------------------------------------
real MeshQuality::mean_ratio(Cell& cell) const
{
  uint const d = cell.dim();
  mapping_.update(cell);

  // Compute the square of the Frobenius norm;
  real sqrFnorm = 0.0;
  for (uint i = 0; i < d; ++i)
  {
    for (uint j = 0; j < d; ++j)
    {
      sqrFnorm += mapping_.J[i + j * Point::MAX_SIZE]
          * mapping_.J[i + j * Point::MAX_SIZE];
    }
  }

  return d * std::pow(std::fabs(mapping_.det), 2.0 / d) / sqrFnorm;
}

//-----------------------------------------------------------------------------
void MeshQuality::compute()
{
  Mesh& m = mesh();

  mu_max = 0.0;
  mu_min = 1.0;

  h_max = 0.0;
  h_min = 1.0e12;

  real mu_sum = 0.0;
  real h_sum = 0.0;

  for (CellIterator c(m); !c.end(); ++c)
  {
    real mu = mean_ratio(*c);
    real h = c->diameter();

    mu_sum += mu;
    h_sum += h;

    mu_max = std::max(mu_max, mu);
    mu_min = std::min(mu_min, mu);

    h_max = std::max(h_max, h);
    h_min = std::min(h_min, h);
  }

  uint const d = m.topology().dim();
  for (uint i = 0; i < d; ++i)
  {
    bbox_min_[i] = 1.0e12;
    bbox_max_[i] = -1.0e12;
  }

  for (VertexIterator vi(m); !vi.end(); ++vi)
  {
    const Vertex& v = *vi;

    Point p = v.point();

    for (uint i = 0; i < d; ++i)
    {
      bbox_min_[i] = std::min(bbox_min_[i], p[i]);
      bbox_max_[i] = std::max(bbox_max_[i], p[i]);
    }
  }

  mu_avg = mu_sum / m.numCells();
  h_avg = h_sum / m.numCells();

  if (m.is_distributed())
  {
    mu_min = reduceMinReal(mu_min);
    mu_max = reduceMaxReal(mu_max);
    mu_avg = reduceAvgReal(mu_avg);

    h_min = reduceMinReal(h_min);
    h_max = reduceMaxReal(h_max);
    h_avg = reduceAvgReal(h_avg);
  }
}
//-----------------------------------------------------------------------------
real MeshQuality::reduceMinReal(real val)
{
  real val_tmp = val;
#ifdef HAVE_MPI
  MPI_Allreduce(&val_tmp, &val, 1, MPI_DOUBLE, MPI_MIN,
                dolfin::MPI::DOLFIN_COMM);
#endif
  return val;
}
//-----------------------------------------------------------------------------
real MeshQuality::reduceMaxReal(real val)
{
  real val_tmp = val;
#ifdef HAVE_MPI
  MPI_Allreduce(&val_tmp, &val, 1, MPI_DOUBLE, MPI_MAX,
                dolfin::MPI::DOLFIN_COMM);
#endif
  return val;
}
//-----------------------------------------------------------------------------
real MeshQuality::reduceAvgReal(real val)
{
  real val_tmp = val;
#ifdef HAVE_MPI
  MPI_Allreduce(&val_tmp, &val, 1, MPI_DOUBLE, MPI_SUM,
                dolfin::MPI::DOLFIN_COMM);
#endif
  return val / dolfin::MPI::numProcesses();
}
//-----------------------------------------------------------------------------
void MeshQuality::disp()
{
  cout << "Mesh quality rank " << dolfin::MPI::processNumber() << ":" << endl;
  cout << "mu_min: " << mu_min << endl;
  cout << "mu_max: " << mu_max << endl;
  cout << "mu_avg: " << mu_avg << endl;
  cout << "h_min: " << h_min << endl;
  cout << "h_max: " << h_max << endl;
  cout << "h_avg: " << h_avg << endl;
}
//-----------------------------------------------------------------------------

}
