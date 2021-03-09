// Copyright (C) 2013-15 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/function/FunctionInterpolation.h>

#include <dolfin/common/maybe_unused.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/fem/ScratchSpace.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/IntersectionDetector.h>
#include <dolfin/mesh/entities/Vertex.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/entities/iterators/VertexIterator.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

void FunctionInterpolation::compute(GenericFunction const& F0, Function& F1)
{
  if (!F1.compatible(F0))
  {
    error("Interpolation between functions with different value shape");
  }
  else if (&F0.mesh() == &F1.mesh())
  {
    interpolateSM(F0, F1);
  }
  else
  {
    interpolateNM(F0, F1);
  }
}

//-----------------------------------------------------------------------------

struct EvalExpression : ufc::function
{
  EvalExpression(Expression const& E) :
    E_(E)
  {
  }

  ///
  inline void evaluate(real* values, const real* coordinates,
                       const ufc::cell&) const override
  {
    E_.eval(values, coordinates);
  }

private:

  Expression const& E_;
};


void FunctionInterpolation::compute(Expression const& F0, Function& F1)
{
  if (!F1.compatible(F0))
  {
    error("FunctionInterpolation: Interpolation between functions with different value shape");
  }
  else
  {
    // Just assume that the coefficient can be evaluated on a ufc::cell
    ScratchSpace S1(F1.space());
    EvalExpression EE(F0);
    size_t dof = 0;
    real * block1 = F1.create_block();
    for (CellIterator cell(F1.mesh()); !cell.end(); ++cell)
    {
      S1.cell.update(*cell);
      // FIXME orientation (0) needs to be correctly set here ?!
      F1.space().element().ufc().evaluate_dofs(&block1[dof], EE,
                                              S1.cell.coordinates.data(),
                                              0, S1.cell);
      dof += S1.local_dimension;
    }
    F1.set_block(block1);
    delete[] block1;
  }
}

//-----------------------------------------------------------------------------
void FunctionInterpolation::compute(Coefficient const& F0, Function& F1)
{
  if (!F1.compatible(F0))
  {
    error("Interpolation between functions with different value shape");
  }
  else
  {
    // Just assume that the coefficient can be evaluated on a ufc::cell
    ScratchSpace S1(F1.space());

    size_t dof = 0;
    real * block1 = F1.create_block();
    for (CellIterator cell(F1.mesh()); !cell.end(); ++cell)
    {
      S1.cell.update(*cell);
      // FIXME orientation (0) needs to be correctly set here ?!
      F1.space().element().ufc().evaluate_dofs(&block1[dof], F0,
                                               S1.cell.coordinates.data(),
                                               0, S1.cell);
      dof += S1.local_dimension;
    }
    F1.set_block(block1);
    delete[] block1;
  }
}

//-----------------------------------------------------------------------------
void FunctionInterpolation::interpolateSM(GenericFunction const& F0,
                                          Function& F1)
{
  message(1, "FunctionInterpolation: interpolation on same mesh");
  dolfin_assert(F0.mesh() == F1.mesh());

  //
  if (F1.space().is_flattenable())
  {
    // Analytical expression and flattened space (naive implementation)
    std::vector<ufc::finite_element const*> const& Sflt =
        F1.space().element().flatten();
    ScratchSpace S1(F1.space());

    size_t dof = 0;
    real * block1 = F1.create_block();
    for (CellIterator cell(F1.mesh()); !cell.end(); ++cell)
    {
      S1.cell.update(*cell);
      S1.finite_element->tabulate_dof_coordinates(S1.coordinates.data(),
                                                  S1.cell.coordinates.data());

      size_t celldof = 0;
      for (size_t leaf = 0; leaf < Sflt.size(); ++leaf)
      {
        for (size_t ii = 0; ii < Sflt[leaf]->space_dimension(); ++ii)
        {
          F0.evaluate(S1.values.data(), &S1.coordinates[Space::MAX_DIMENSION*(celldof++)], S1.cell);
          block1[dof++] = S1.values[leaf];
        }
      }
      //
      dolfin_assert(celldof == S1.local_dimension);
    }
    F1.set_block(block1);
    delete[] block1;
  }
  else
  {
    // The other function is discrete and non-trivial
    ScratchSpace S1(F1.space());

    size_t dof = 0;
    real * block1 = F1.create_block();
    for (CellIterator cell(F1.mesh()); !cell.end(); ++cell)
    {
      S1.cell.update(*cell);
      // FIXME orientation (0) needs to be correctly set here ?!
      F1.space().element().ufc().evaluate_dofs(&block1[dof], F0,
                                               S1.cell.coordinates.data(),
                                               0, S1.cell);
      dof += S1.local_dimension;
    }
    F1.set_block(block1);
    delete[] block1;
  }
}

//-----------------------------------------------------------------------------
void FunctionInterpolation::interpolateNM(GenericFunction const& F0,
                                          Function& F1)
{
  message(1, "FunctionInterpolation: interpolation on non-matching meshes");
  dolfin_assert(F0.mesh() != F1.mesh());

  Mesh& M0 = F0.mesh();
  //size_t const gdim0 = M0.geometry_dimension();
  //size_t const tdim0 = M0.geometry_dimension();
  Cell c00(M0,0);
  UFCCell ufc0(c00);

  Mesh& M1 = F1.mesh();
  size_t const gdim1 = M1.geometry_dimension();
  size_t const tdim1 = M1.geometry_dimension();
  FiniteElementSpace const& Vh1 = F1.space();
  DofMap const& dm1 = Vh1.dofmap();
  ScratchSpace S1(Vh1);

  //
#if DOLFIN_HAVE_MPI
  size_t rank = dolfin::MPI::rank();
#endif
  size_t pe_size = dolfin::MPI::size();

  // On-proc for M0 and M1
  std::vector<size_t> dofs_indices0;
  std::vector<size_t> dofs_valsidx0;
  std::vector<size_t> cell_indices0;
  std::vector<real> dofs_xcoords0;

  // Off-proc for M0: need to recv
  std::vector<size_t> dofs_indicesX;
  std::vector<real> dofs_xcoordsX;

  // Off-proc for M1: need to send
  std::vector<size_t> dofs_indices1;
  std::vector<size_t> dofs_valsidx1;
  std::vector<size_t> cell_indices1;
  std::vector<real> dofs_xcoords1;

  // DEBUG
  _set<size_t> offproc;

  // Dofs count to be sent and received
#if DOLFIN_HAVE_MPI
  size_t num_sendadj = 0;
  size_t num_recvadj = 0;
#endif

  std::vector< size_t > dof1sendcount( pe_size, 0 );
  std::vector< size_t > dof1recvcount( pe_size, 0 );

  // Total count to be received from other ranks
  size_t num_dofsF = 0;

  // Some flags
  bool const is_distributed = M0.is_distributed() || M1.is_distributed();
#if DOLFIN_HAVE_MPI
  bool const just_first_coords = Vh1.is_flattenable()
                                 && Vh1.element().is_vectorizable();
#endif

  //--- Collect on-proc and off-proc dofs
  if (Vh1.is_vertex_based())
  {
    // This implementation assumes a scalar or vector-valued function which has
    // dofs only located at vertices.
    // (u, r) : (dof indices located at vertex, vertex coordinates)
    size_t const num_cellverts = M1.type().num_entities(0);
    M1.init(0, tdim1);
    for (VertexIterator v1(M1); !v1.end(); ++v1)
    {
      if (!v1->is_ghost())
      {
        Cell c1(M1, v1->entities(tdim1)[0]);
        S1.cell.update(c1);
        dm1.tabulate_dofs(S1.dofs.data(), S1.cell);
        std::vector<size_t> & vid = c1.entities(0);
        size_t vpos = 0;
        while (vid[vpos] != v1->index())
        {
          ++vpos;
        }
        Point p = v1->point();
        std::vector<size_t> M0cells;
        M0.intersector().overlap(p, M0cells);
        if (M0cells.empty())
        {
          // Global dof indices
          for (size_t v = 0; v < S1.size; ++v)
          {
            dofs_indicesX.push_back(S1.dofs[vpos + v * num_cellverts]);
          }
          // Coordinates
          for (size_t d = 0; d < gdim1; ++d)
          {
            dofs_xcoordsX.push_back(p[d]);
          }
        }
        else
        {
          // Check consistency
          dolfin_assert(M0.type().intersects(Cell(M0, M0cells.front()), p));
          // Local cell
          cell_indices0.push_back(M0cells.front());
          // Global dof indices
          for (size_t v = 0; v < S1.size; ++v)
          {
            dofs_indices0.push_back(S1.dofs[vpos + v * num_cellverts]);
          }
          // Coordinates
          for (size_t d = 0; d < gdim1; ++d)
          {
            dofs_xcoords0.push_back(p[d]);
          }
        }
      }
    }
    // DEBUG
    dolfin_assert(dofs_indicesX.size() / S1.size == dofs_xcoordsX.size() / gdim1);
    //
    dolfin_assert(dofs_indices0.size() / S1.size == dofs_xcoords0.size() / gdim1);
  }
  else if (Vh1.is_cellwise_constant())
  {
    // This implementation assumes a scalar or vector-valued piecewise constant
    // function which has dofs only located at the cell midpoint
    // (u, r) : (dof indices located at cell c, cell midpoint coordinates)
    for (CellIterator c1(M1); !c1.end(); ++c1)
    {
      S1.cell.update(*c1);
      dm1.tabulate_dofs(S1.dofs.data(), S1.cell);
      Point p = c1->midpoint();
      std::vector<size_t> M0cells;
      M0.intersector().overlap(p, M0cells);
      if (M0cells.empty())
      {
        // Global dof indices
        for (size_t v = 0; v < S1.size; ++v)
        {
          dofs_indicesX.push_back(S1.dofs[v]);
        }
        // Coordinates
        for (size_t d = 0; d < gdim1; ++d)
        {
          dofs_xcoordsX.push_back(p[d]);
        }
      }
      else
      {
        // Check consistency
        dolfin_assert(M0.type().intersects(Cell(M0, M0cells.front()), p));
        // Local cell
        cell_indices0.push_back(M0cells.front());
        // Global dof indices
        for (size_t v = 0; v < S1.size; ++v)
        {
          dofs_indices0.push_back(S1.dofs[v]);
        }
        // Coordinates
        for (size_t d = 0; d < gdim1; ++d)
        {
          dofs_xcoords0.push_back(p[d]);
        }
      }
    }
    // DEBUG
    dolfin_assert(dofs_indicesX.size() / S1.size == dofs_xcoordsX.size() / gdim1);
    //
    dolfin_assert(dofs_indices0.size() / S1.size == dofs_xcoords0.size() / gdim1);
  }
  else if (Vh1.is_flattenable() && Vh1.element().is_vectorizable())
  {
    // This implementation assumes a scalar function for which
    // components are approximated in a discrete space other than CG1 and DG0.
    // (u, r) : (indices for dofs located at node n,
    //           node coordinates)
    _set<size_t> done;
    size_t const local_dim1 = dm1.num_element_dofs / S1.size;
    Point p;
    for (CellIterator c1(M1); !c1.end(); ++c1)
    {
      S1.cell.update(*c1);
      dm1.tabulate_dofs(S1.dofs.data(), S1.cell);
      S1.finite_element->tabulate_dof_coordinates(S1.coordinates.data(),
                                                  S1.cell.coordinates.data());

      // For each dof of the first leaf
      for (size_t i = 0; i < local_dim1; ++i)
      {
        if ((done.count(S1.dofs[i]) == 0) && !dm1.is_ghost(S1.dofs[i]))
        {
          done.insert(S1.dofs[i]);
          std::memcpy(&p[0], &S1.coordinates[i*Space::MAX_DIMENSION], gdim1 * sizeof(real));
          std::vector<size_t> M0cells;
          M0.intersector().overlap(p, M0cells);
          if (M0cells.empty())
          {
            // Global dof indices
            for (size_t v = 0; v < S1.size; ++v)
            {
              dofs_indicesX.push_back(S1.dofs[v * local_dim1 + i]);
            }
            // Coordinates
            for (size_t d = 0; d < gdim1; ++d)
            {
              dofs_xcoordsX.push_back(p[d]);
            }
          }
          else
          {
            // Check consistency
            dolfin_assert(M0.type().intersects(Cell(M0,M0cells.front()), p));
            // Local cell
            cell_indices0.push_back(M0cells.front());
            // Global dof indices
            for (size_t v = 0; v < S1.size; ++v)
            {
              dofs_indices0.push_back(S1.dofs[v * local_dim1 + i]);
            }
            // coordinates
            for (size_t d = 0; d < gdim1; ++d)
            {
              dofs_xcoords0.push_back(p[d]);
            }
          }
        }
      }
    }
  }
  else
  {
    error("Unsupported discrete space");
  }

  // DEBUG
  for (size_t i = 0; i < dofs_indicesX.size(); ++i)
  {
    offproc.insert(dofs_indicesX[i]);
  }

  //--- Exchange off-proc dofs
  if (is_distributed)
  {
#ifdef DOLFIN_HAVE_MPI

    int u_maxrecvcount = dofs_indicesX.size();
    MPI::all_reduce_in_place<MPI::max>( u_maxrecvcount );
    int r_maxrecvcount = u_maxrecvcount * gdim1;

    //
    std::vector< size_t > u_recvbuf( u_maxrecvcount );
    std::vector< real > r_recvbuf( r_maxrecvcount );
    for (int j = 1; j < (int) pe_size; ++j)
    {
      int src = (rank - j + pe_size) % pe_size;
      int dest = (rank + j) % pe_size;

      int u_recvcount = MPI::sendrecv( dofs_indicesX, dest, u_recvbuf, src, 1 );
      int r_recvcount = MPI::sendrecv( dofs_xcoordsX, dest, r_recvbuf, src, 1 );

      size_t matching_dofs = 0;
      if (u_recvcount > 0)
      {
        if (just_first_coords)
        {
          // DEBUG
          dolfin_assert(u_recvcount / S1.size == r_recvcount / gdim1);
          MAYBE_UNUSED(r_recvcount);

          //
          Point p;
          size_t node_count = u_recvcount / S1.size;
          for (size_t i = 0; i < node_count; ++i)
          {
            std::vector<size_t> M0cells;
            std::memcpy(&p[0], &r_recvbuf[i * gdim1], gdim1 * sizeof(real));
            M0.intersector().overlap(p, M0cells);
            if (!M0cells.empty())
            {
              // Check consistency
              dolfin_assert(M0.type().intersects(Cell(M0, M0cells.front()), p));
              // Local cell
              cell_indices1.push_back(M0cells.front());
              // Global dof indices
              for (size_t v = 0; v < S1.size; ++v)
              {
                dofs_indices1.push_back(u_recvbuf[i * S1.size + v]);
                ++matching_dofs;
              }
              // Coordinates
              for (size_t d = 0; d < gdim1; ++d)
              {
                dofs_xcoords1.push_back(p[d]);
              }
            }
          }
          //
          dolfin_assert(dofs_indices1.size() == cell_indices1.size() * S1.size);
          //
          dolfin_assert(cell_indices1.size() * gdim1 == dofs_xcoords1.size());
        }
        else
        {
          error("Only scalar and vector valued functions are supported.");
        }
      }
      // Set number of dofs to be sent back to src
      dof1sendcount[src] = matching_dofs;
    }

    // Collect data on dofs distribution across ranks
    for (int j = 1; j < (int) pe_size; ++j)
    {
      int dest = (rank + j) % pe_size;
      int src = (rank - j + pe_size) % pe_size;

      MPI::sendrecv(&dof1sendcount[dest], 1, dest,
                    &dof1recvcount[src], 1, src, 1 );

      //
      if (dof1sendcount[dest] > 0)
      {
        ++num_sendadj;
      }
      if (dof1recvcount[src] > 0)
      {
        ++num_recvadj;
      }
      num_dofsF += dof1recvcount[src];

    }

#if DEBUG
    message("Rank %4d: In send order", rank);
    for (int j = 1; j < (int) pe_size; ++j)
    {
      int dest = (rank - j + pe_size) % pe_size;
      message("Rank %4d: Send = %8d Recv = %8d", dest, dof1sendcount[dest],
              dof1recvcount[dest]);
    }
    message("TOTAL    : Send = %8d Recv = %8d", dofs_indices1.size(),
            num_dofsF);
#endif
    dolfin_assert(dof1sendcount[rank] == 0);
#endif
  }

  //--- Evaluation

  // Prepare requests to receive dof indices and values
  std::vector< size_t > dofs_indicesF( num_dofsF );
  std::vector< real > dofs_cvaluesF( num_dofsF );

#ifdef DOLFIN_HAVE_MPI
  MPI_Status status;
  std::vector< MPI_Request > u_req_recv( num_recvadj );
  std::vector< MPI_Request > r_req_recv( num_recvadj );
  if (is_distributed)
  {
    size_t offsetF = 0;
    size_t recv_id = 0;
    for (int j = 1; j < (int) pe_size; ++j)
    {
      int src = (rank + j) % pe_size;
      size_t count = dof1recvcount[src];
      if (count > 0)
      {
        MPI::check_error( MPI_Irecv(&dofs_indicesF[offsetF], count,
                                    MPI_type< size_t >::value, src, 0,
                                    MPI::DOLFIN_COMM, &u_req_recv[recv_id]) );
        MPI::check_error( MPI_Irecv(&dofs_cvaluesF[offsetF], count,
                                    MPI_type< real >::value, src, 0,
                                    MPI::DOLFIN_COMM, &r_req_recv[recv_id]) );
        offsetF += count;
        ++recv_id;
      }
    }
  }
#endif

  // Prepare requests to send dof indices and values
  std::vector< real > dofs_cvalues1( dofs_indices1.size() );
#ifdef DOLFIN_HAVE_MPI
  std::vector< MPI_Request > u_req_send( num_sendadj );
  std::vector< MPI_Request > r_req_send( num_sendadj );
  if (is_distributed)
  {
    size_t offset1 = 0;
    size_t send_id = 0;
    if (just_first_coords)
    {
      for (int j = 1; j < (int) pe_size; ++j)
      {
        int dest = (rank - j + pe_size) % pe_size;
        size_t count = dof1sendcount[dest];
        if (count > 0)
        {
          size_t cell_offset = offset1 / S1.size;
          size_t node_count = count / S1.size;
          Point p;
          for (size_t nodei = cell_offset; nodei < (cell_offset + node_count);
              ++nodei)
          {
            Cell c0(M0, cell_indices1[nodei]);
            ufc0.update(c0);

            // Check consistency
            std::memcpy(&p[0], &dofs_xcoords1[nodei * gdim1],
                        sizeof(real) * gdim1);
            dolfin_assert(M0.type().intersects(c0, p));

            // Evaluate
            F0.evaluate(&dofs_cvalues1[nodei * S1.size], &p[0], ufc0);
          }
          MPI::check_error( MPI_Isend(&dofs_indices1[offset1], count,
                                      MPI_type< size_t >::value, dest, 0,
                                      MPI::DOLFIN_COMM, &u_req_send[send_id]) );
          MPI::check_error( MPI_Isend(&dofs_cvalues1[offset1], count,
                                      MPI_type< real >::value, dest, 0,
                                      MPI::DOLFIN_COMM, &r_req_send[send_id]) );
          offset1 += count;
          ++send_id;
        }
      }
    }
  }
#endif

  // Local dofs
  Point n;
  size_t const num_dofs0 = dofs_indices0.size();
  size_t const num_node0 = cell_indices0.size();
  real * dofs_cvalues0 = new real[num_dofs0];
  for (size_t ii = 0; ii < num_node0; ++ii)
  {
    Cell c0(M0, cell_indices0[ii]);
    ufc0.update(c0);

    // Check consistency
    std::memcpy(&n[0], &dofs_xcoords0[ii * gdim1], sizeof(real) * gdim1);
    dolfin_assert(M0.type().intersects(c0, n));

    // Let us just use a point for now
    F0.evaluate(&dofs_cvalues0[ii * S1.size], &n[0], ufc0);
  }
  F1.vector().set(&dofs_cvalues0[0], dofs_indices0.size(), &dofs_indices0[0]);
  delete[] dofs_cvalues0;

  // End off proc
  if (is_distributed)
  {
#ifdef DOLFIN_HAVE_MPI
    for (int j = 0; j < (int) num_sendadj; ++j)
    {
      MPI::check_error( MPI_Wait(&u_req_send[j], &status) );
      MPI::check_error( MPI_Wait(&r_req_send[j], &status) );
    }
    for (int j = 0; j < (int) num_recvadj; ++j)
    {
      MPI::check_error( MPI_Wait(&u_req_recv[j], &status) );
      MPI::check_error( MPI_Wait(&r_req_recv[j], &status) );
    }
#endif
  }

  // DEBUG
  for (size_t i = 0; i < num_dofsF; ++i)
  {
    if (offproc.count(dofs_indicesF[i]) == 0)
    {
      error("Trying to set invalid dof %d", dofs_indicesF[i]);
    }
  }

  // Set foreign dofs values
  F1.vector().set(&dofs_cvaluesF[0], num_dofsF, &dofs_indicesF[0]);
  F1.sync();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
