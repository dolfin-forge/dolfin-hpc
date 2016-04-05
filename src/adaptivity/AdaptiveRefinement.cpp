// Copyright (C) 2010 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-09-13
// Last changed: 2011-04-16

#include <algorithm>
#include <limits>
#include <fstream>
#include <errno.h>
#include <sys/stat.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/io/BinaryFile.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionDecomposition.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/LoadBalancer.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/RivaraRefinement.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/parameter/parameters.h>

#include <dolfin/adaptivity/AdaptiveRefinement.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace dolfin
{

//-----------------------------------------------------------------------------
void AdaptiveRefinement::refine(Mesh& mesh, MeshFunction<bool>& cell_marker)

{
  message("Adaptive refinement");

  uint const numcellsbefore = mesh.global_numCells();
  uint const numvertsbefore = mesh.global_numVertices();
  message("  - cells    before: %d", numcellsbefore);
  message("  - vertices before: %d", numvertsbefore);

  std::ostringstream marked_filename;
  marked_filename << "marked";
  std::string const marked_format = dolfin_get("output_format");
  if (marked_format == "vtk")
  {
    marked_filename << ".pvd";
  }
  else if (marked_format == "binary")
  {
    marked_filename << ".bin";
  }
  File refinefile(marked_filename.str());
  refinefile << cell_marker;

  std::string const refine_type = dolfin_get("adapt_algorithm");
  if (refine_type == "rivara")
  {
    RivaraRefinement::refine(mesh, cell_marker);
  }
  else if (refine_type == "simple")
  {
    mesh.refine(cell_marker, true);
  }
  else
  {
    dolfin::error("Unknown refinement algorithm");
  }

  uint const numcellsafter = mesh.global_numCells();
  dolfin_assert(numcellsafter >= numcellsbefore);
  uint const numvertsafter = mesh.global_numVertices();
  dolfin_assert(numvertsafter >= numvertsbefore);
  message("  - cells    after: %d", numcellsafter);
  message("  - vertices after: %d", numvertsafter);
  skip();
}
//-----------------------------------------------------------------------------
void AdaptiveRefinement::refine_and_project(Mesh& mesh,
                                            Array<Function *> const& functions,
                                            MeshFunction<bool>& cell_marker)
{

  dolfin_set("Load balancer redistribute", false);
  message("Adaptive refinement (with projection)");
  message("  - cells    before: %d", mesh.global_numCells());
  message("  - vertices before: %d", mesh.global_numVertices());

  std::string const refine_type = dolfin_get("adapt_algorithm");
  if (refine_type == "simple")
  {
    LoadBalancer::balance(mesh, cell_marker);
  }
  else if (refine_type == "rivara")
  {
    LoadBalancer::balance(mesh, cell_marker, LoadBalancer::LEPP);
  }
  else
  {
    dolfin::error("Unknown refinement algorithm");
  }

  std::ostringstream marked_filename;
  marked_filename << "marked";
  std::string const marked_format = dolfin_get("output_format");
  if (marked_format == "vtk")
  {
    marked_filename << ".pvd";
  }
  else if (marked_format == "binary")
  {
    marked_filename << ".bin";
  }
  File refinefile(marked_filename.str());
  refinefile << cell_marker;

  error("Fix the shit");
  MeshFunction<uint> *partitions = NULL;// mesh.data().meshFunction("partitions");

  uint const maxvecsize = 3;
  real * x_values[maxvecsize];
  uint * x_rows[maxvecsize];
  uint x_m[maxvecsize];

  for (Array<Function *>::const_iterator it = functions.begin();
      it != functions.end(); ++it)
  {
    Function const& func_to_project = **it;
    if (func_to_project.type() != Function::discrete)
    {
      error("Projection only implemented for discrete functions");
    }

    FiniteElementSpace const& space = func_to_project.space();
    uint const num_sub = space.element().num_sub_elements();

    if (num_sub == 0)
    {
      continue;
    }

    Array<Function *> coarse = FunctionDecomposition::compute(func_to_project);

    for (uint i = 0; i < num_sub; ++i)
    {
      AdaptiveRefinement::redistribute_func(mesh, *coarse[i], &x_values[i],
                                            &x_rows[i], x_m[i], *partitions);
    }

    while (!coarse.empty())
    {
      delete coarse.back();
      coarse.pop_back();
    }
  }

  MeshFunction<bool> new_cell_marker;
  mesh.distribute(*partitions, cell_marker, new_cell_marker);

  Mesh new_mesh = mesh;
  RivaraRefinement::refine(new_mesh, new_cell_marker, 0.0, 0.0, 0.0, false);
  new_mesh.renumber();

  if (MPI::processNumber() == 0)
  {
    if (mkdir("../scratch", S_IRWXU) < 0)
    {
      perror("mkdir failed");
    }
  }

  int p_count = 0;
  for (Array<Function *>::const_iterator it = functions.begin();
      it != functions.end(); ++it)
  {
    Function const& func_to_project = **it;
    if (func_to_project.type() != Function::discrete)
    {
      error("Projection only implemented for discrete functions");
    }

    Array<Function *> post;
    FiniteElementSpace const& space = func_to_project.space();
    uint const num_sub = space.element().num_sub_elements();

    // FIXME: Invalid for scalar functions due to the zero subspace assumption
    for (uint i = 0; i < num_sub; ++i)
    {
      FiniteElementSpace subspace(space, i);
      post.push_back(new Function(subspace));
      post.back()->vector().set(x_values[i], x_m[i], x_rows[i]);
      post.back()->sync_ghosts();
    }

    for (uint i = 0; i < num_sub; ++i)
    {
      delete[] x_values[i];
      delete[] x_rows[i];
    }

    FiniteElementSpace projected_space(new_mesh, space);
    Function proj(projected_space);
    AdaptiveRefinement::project(new_mesh, post, proj);

    while (!post.empty())
    {
      delete post.back();
      post.pop_back();
    }

    std::stringstream p_filename;
#ifdef ENABLE_MPIIO
    p_filename << "../scratch/projected_" << p_count++ << ".bin" << std::ends;
#else
    p_filename << "../scratch/projected_" << p_count++ << "_" << MPI::processNumber() << ".bin" << std::ends;
#endif
    File p_file(p_filename.str());
    p_file << proj.vector();
  }

  mesh = new_mesh;
  mesh.renumber();

}
//-----------------------------------------------------------------------------
void AdaptiveRefinement::redistribute_func(Mesh& mesh, Function const& f,
                                           real **vp, uint **rp, uint& m,
                                           MeshFunction<uint>& distribution)
{

#ifdef HAVE_MPI
  uint pe_rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint target_proc, src, dest, recv_size, local_size;

  MPI_Status status;

  real *values = new real[f.vector().local_size()];
  f.vector().get(values);

  Array<real> *send_buffer = new Array<real> [pe_size];
  Array<uint> *send_buffer_indices = new Array<uint> [pe_size];

  MeshFunction<bool> marked(mesh, 0);
  marked = false;

  std::vector<std::pair<uint, real> > recv_data;
  local_size = 0;

  real value;
  uint global_index;

  for (CellIterator c(mesh); !c.end(); ++c)
  {

    target_proc = distribution.get(*c);

    //FIXME: Only P1 friendly.
    for (VertexIterator v(*c); !v.end(); ++v)
    {

      global_index = mesh.distdata()[0].get_global(v->index());
      f.vector().get(&value, 1, &global_index);

      if (target_proc == pe_rank && !v->is_ghost()
          && !marked.get(*v))
      {

        std::pair<uint, real> p(global_index, value);
        recv_data.push_back(p);
        marked.set(*v, true);
        continue;
      }

      if (!v->is_ghost() && !marked.get(*v))
      {

        send_buffer[target_proc].push_back(value);
        send_buffer_indices[target_proc].push_back(global_index);
        marked.set(*v, true);
      }
    }

    local_size = std::max(local_size, (uint) send_buffer[target_proc].size());
  }

  delete[] values;

  MPI_Allreduce(&local_size, &recv_size, 1, MPI_UNSIGNED, MPI_MAX,
                MPI::DOLFIN_COMM);

  real *recv_buffer = new real[recv_size];
  uint *recv_buffer_indices = new uint[recv_size];

  //
  int recv_count = 0;
  for (uint j = 1; j < pe_size; j++)
  {
    src = (pe_rank - j + pe_size) % pe_size;
    dest = (pe_rank + j) % pe_size;

    MPI_Sendrecv(&send_buffer[dest][0], send_buffer[dest].size(), MPI_DOUBLE,
                 dest, 1, recv_buffer, recv_size, MPI_DOUBLE, src, 1,
                 MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&send_buffer_indices[dest][0],
                 send_buffer_indices[dest].size(), MPI_UNSIGNED, dest, 1,
                 recv_buffer_indices, recv_size, MPI_UNSIGNED, src, 1,
                 MPI::DOLFIN_COMM, &status);

    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    for (int i = 0; i < recv_count; i++)
    {
      std::pair<uint, real> p(recv_buffer_indices[i], recv_buffer[i]);
      recv_data.push_back(p);
    }

  }

  delete[] recv_buffer_indices;
  delete[] recv_buffer;

  for (uint i = 0; i < pe_size; i++)
  {
    send_buffer[i].clear();
    send_buffer_indices[i].clear();
  }

  delete[] send_buffer_indices;
  delete[] send_buffer;

  local_size = recv_data.size();
  values = new real[local_size];
  uint *rows = new uint[local_size];

  for (uint i = 0; i < local_size; i++)
  {
    rows[i] = recv_data[i].first;
    values[i] = recv_data[i].second;
  }

  *vp = values;
  *rp = rows;
  m = local_size;

#endif

}
//-----------------------------------------------------------------------------
void AdaptiveRefinement::project(Mesh& new_mesh, Array<Function *>& f_post,
                                 Function& projected)
{
  GenericVector& x_proj = projected.vector();
  real test_value;
  real x[3];
  real *vv = new real[x_proj.local_size()];
  uint *indices = new uint[x_proj.local_size()];
  uint *local_indices = new uint[projected.space().dofmap().local_dimension()];
  uint i = 0;
  MeshFunction<bool> processed(new_mesh, 0);
  processed = false;

  projected.vector().zero();
  projected.sync_ghosts();

  real gts_tol = dolfin_get("GTS Tolerance");
  real geom_tol = dolfin_get("Geometrical Tolerance Tetrahedron");

  dolfin_set("GTS Tolerance", 1e-10);
  dolfin_set("Geometrical Tolerance Tetrahedron", 1e-8);
  CellIterator c(new_mesh);
  UFCCell ufccell(*c);
  if ((projected.space().family() == ufl::Family::CG)
      && (projected.space().degree() == 1))
  {
    for (; !c.end(); ++c)
    {

      ufccell.update(*c);
      projected.space().dofmap().tabulate_dofs(local_indices, ufccell, c->index());

      //FIXME: Only P1 friendly.
      for (VertexIterator v(*c); !v.end(); ++v)
      {

        uint *cvi = c->entities(0);
        uint ci = 0;
        for (ci = 0; ci < c->num_entities(0); ci++)
        {
          if (cvi[ci] == v->index())
          {
            break;
          }
        }

        if (v->is_ghost() || processed.get(*v))
        {
          continue;
        }
        processed.set(*v, true);

        Vertex *v_e = 0;

        x[0] = v->x()[0];
        x[1] = v->x()[1];
        x[2] = v->x()[2];
        f_post[0]->eval(&test_value, &x[0]);
        if (test_value == std::numeric_limits<real>::infinity())
        {
          for (EdgeIterator e(*v); !e.end(); ++e)
          {
            uint const *edge_v = e->entities(0);

            if (edge_v[0] != v->index())
            {
              v_e = new Vertex(new_mesh, edge_v[0]);
            }
            else
            {
              v_e = new Vertex(new_mesh, edge_v[1]);
            }

            x[0] = v_e->x()[0];
            x[1] = v_e->x()[1];
            x[2] = v_e->x()[2];
            f_post[0]->eval(&test_value, &x[0]);

            if (test_value != std::numeric_limits<real>::infinity())
            {
              break;
            }

          }

          if (test_value == std::numeric_limits<real>::infinity())
          {
            error("Couldn't find any suitable projection point");
          }
        }

        vv[i] = test_value;
        indices[i++] = local_indices[ci];

        for (uint d = 1; d < f_post.size(); ++d)
        {
          f_post[d]->eval(&test_value, &x[0]);
          if (test_value != std::numeric_limits<real>::infinity())
          {
            vv[i] = test_value;
            indices[i++] = local_indices[ci + d * c->num_entities(0)];
          }
        }

        delete v_e;
      }
    }
  }
  else
  {
    error("AdaptiveRefinement::project only implemented for P1");
  }
  dolfin_set("GTS Tolerance", gts_tol);
  dolfin_set("Geometrical Tolerance Tetrahedron", geom_tol);

  if (i > 0)
  {
    x_proj.set(vv, i, indices);
    x_proj.apply();
  }
  projected.sync_ghosts();

  delete[] vv;
  delete[] indices;
  delete[] local_indices;
}
//-----------------------------------------------------------------------------

}

