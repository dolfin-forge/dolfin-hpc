// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Magnus Vikstrøm, 2007.
// Modified by Anders Logg, 2007.
// Modified by Niclas Jansson, 2008-2015.
// Modified by Balthasar Reuter, 2013.
//
// First added:  2007-05-30
// Last changed: 2015-01-31

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MPIMeshCommunicator.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <map>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace dolfin;

#ifdef HAVE_MPI

//-----------------------------------------------------------------------------
MPIMeshCommunicator::MPIMeshCommunicator()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MPIMeshCommunicator::~MPIMeshCommunicator()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(Mesh& mesh,
                                     MeshFunction<uint>& distribution)
{
  distributeCommon(mesh, distribution, (MeshFunction<bool> *) 0,
                   (MeshFunction<bool> *) 0);
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(Mesh& mesh,
                                     MeshFunction<uint>& distribution,
                                     MeshFunction<bool>& old_cell_marker,
                                     MeshFunction<bool>& cell_marker)
{
  distributeCommon(mesh, distribution, &old_cell_marker, &cell_marker);
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(
    Mesh& mesh,
    MeshFunction<uint>& distribution,
    Array<std::pair<MeshFunction<uint> *, MeshFunction<uint> *> >& cell_functions)
{
  distributeCommon(mesh, distribution, &cell_functions, 0);
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(
    Mesh& mesh,
    MeshFunction<uint>& distribution,
    Array<std::pair<MeshFunction<real> *, MeshFunction<real> *> >& vertex_functions)
{
  distributeCommon(mesh, distribution, 0, &vertex_functions);
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(
    Mesh& mesh,
    MeshFunction<uint>& distribution,
    Array<std::pair<MeshFunction<uint> *, MeshFunction<uint> *> >& cell_functions,
    Array<std::pair<MeshFunction<real> *, MeshFunction<real> *> >& vertex_functions)
{
  distributeCommon(mesh, distribution, &cell_functions, &vertex_functions);
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distributeCommon(Mesh& mesh,
                                           MeshFunction<uint>& distribution,
                                           MeshFunction<bool> *old_cell_marker,
                                           MeshFunction<bool> *cell_marker)
{

  MeshDistributedData distdata(mesh.topology().dim());
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint gdim = mesh.geometry().dim();
  uint ndims = mesh.type().num_vertices(mesh.topology().dim());

  Array<real> *send_list_vertices = new Array<real> [pe_size];
  Array<uint> *send_list_mappings = new Array<uint> [pe_size];
  Array<uint> *send_list_cells = new Array<uint> [pe_size];
  Array<real> coords;
  Array<uint> cl, shared_buffer;
  Array<bool> cm;
  uint num_cells, num_vertices, target_proc, offset;

  int recv_size, recv_size_cell, send_size;
  int recv_count, recv_count_vertices, recv_count_cells;

  // Process mesh entities according to distribution
  uint vi = 0;
  // Distribution defined per vertex
  if (distribution.dim() == 0)
  {
    for (VertexIterator v(mesh); !v.end(); ++v)
    {
      if (!v->is_ghost())
      {
        if (distribution.get(*v) != rank)
        {
          target_proc = distribution.get(*v);
          send_list_mappings[target_proc].push_back(v->global_index());
          for (uint d = 0; d < gdim; ++d)
          {
            send_list_vertices[target_proc].push_back(v->x(d));
          }
        }
        else
        {
          for (uint d = 0; d < gdim; ++d)
          {
            coords.push_back(v->x(d));
          }
          distdata[0].set_map(vi++, v->global_index());
        }
      }
    }
    recv_count_cells = 0;
  }
  // Distribution defined per cell
  else if (distribution.dim() == mesh.topology().dim())
  {

    MeshFunction<bool> vertex_used(mesh, 0);
    vertex_used = false;

    for (CellIterator c(mesh); !c.end(); ++c)
    {
      if (distribution.get(*c) != rank)
      {
        target_proc = distribution.get(*c);
        for (VertexIterator v(*c); !v.end(); ++v)
        {
          // Buffer cell global vertex indices
          uint const glb_index = v->global_index();
          send_list_cells[target_proc].push_back(glb_index);

          // Buffer all cell vertices that belong to another processor
          if (!v->is_ghost() && !vertex_used.get(*v))
          {
            send_list_mappings[target_proc].push_back(glb_index);
            for (uint d = 0; d < gdim; ++d)
            {
              send_list_vertices[target_proc].push_back(v->x(d));
            }
            vertex_used.set(*v, true);
          }
        }

        // Transfer Cell marker  (mesh refinement)
        if (cell_marker)
        {
          if (old_cell_marker->get(*c))
          {
            send_list_cells[target_proc].push_back(1);
          }
          else
          {
            send_list_cells[target_proc].push_back(0);
          }
        }
      }
      else
      {
        for (VertexIterator v(*c); !v.end(); ++v)
        {
          if (!vertex_used.get(*v))
          {
            if (!v->is_ghost())
            {
              for (uint d = 0; d < gdim; ++d)
              {
                coords.push_back(v->x(d));
              }
              distdata[0].set_map(vi++, v->global_index());
              vertex_used.set(*v, true);
            }
          }
        }
      }
    }

    recv_count_cells = 0;
    for (uint i = 0; i < pe_size; i++)
    {
      send_size = send_list_cells[i].size();
      MPI_Reduce(&send_size, &recv_count_cells, 1, MPI_INT, MPI_SUM, i,
                 MPI::DOLFIN_COMM);
    }
  }
  else
  {
    error("Distribution defined on unkown mesh entity");
  }

  // Exchange the processed entities
  recv_count = 0;
  for (uint i = 0; i < pe_size; i++)
  {
    send_size = send_list_vertices[i].size();
    MPI_Reduce(&send_size, &recv_count, 1, MPI_INT, MPI_MAX, i,
               MPI::DOLFIN_COMM);
  }
  recv_count_vertices = recv_count / gdim;
  num_vertices = recv_count;

  double *recv_buff = new double[recv_count];
  uint *recv_buff_map = new uint[recv_count_vertices];

  num_cells = recv_count_cells;
  uint *recv_buff_cell = new uint[recv_count_cells];
  uint *rcp = &recv_buff_cell[0];

  MPI_Status status;
  uint src, dest, buff_map;
  for (uint i = 1; i < pe_size; i++)
  {

    src = (rank - i + pe_size) % pe_size;
    dest = (rank + i) % pe_size;

    MPI_Sendrecv(&send_list_vertices[dest][0], send_list_vertices[dest].size(),
                 MPI_DOUBLE, dest, 0, recv_buff, recv_count, MPI_DOUBLE, src, 0,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &recv_size);

    MPI_Sendrecv(&send_list_mappings[dest][0], send_list_mappings[dest].size(),
                 MPI_UNSIGNED, dest, 1, recv_buff_map, recv_count_vertices,
                 MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&send_list_cells[dest][0], send_list_cells[dest].size(),
                 MPI_UNSIGNED, dest, 2, rcp, recv_count_cells, MPI_UNSIGNED,
                 src, 2, MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_size_cell);
    rcp += recv_size_cell;
    recv_count_cells -= recv_size_cell;

    buff_map = 0;
    for (int i = 0; i < recv_size; i += gdim)
    {
      if (!distdata[0].has_global(recv_buff_map[buff_map]))
      {
        distdata[0].set_map(vi++, recv_buff_map[buff_map]);
        for (uint d = 0; d < gdim; ++d)
        {
          coords.push_back(recv_buff[i + d]);
        }
      }
      buff_map++;
    }
  }

  //Clear send buffers
  for (uint i = 0; i < pe_size; i++)
  {
    send_list_cells[i].clear();
    send_list_vertices[i].clear();
    send_list_mappings[i].clear();
  }
  delete[] send_list_vertices;
  delete[] send_list_mappings;
  delete[] send_list_cells;
  delete[] recv_buff_map;
  delete[] recv_buff;

  // Process new and old cells if distribution is defined on cells
  if (distribution.dim() == mesh.topology().dim())
  {

    // Add old cells
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      if (distribution.get(*c) == rank)
      {
        for (VertexIterator v(*c); !v.end(); ++v)
        {
          uint const glb_index = v->global_index();
          if (!distdata[0].has_global(glb_index))
          {
            cl.push_back(vi);
            for (uint j = 0; j < gdim; ++j)
            {
              coords.push_back(0.0);
            }
            distdata[0].set_map(vi, glb_index);
            distdata[0].set_ghost(vi++, distribution.get(*c));
            shared_buffer.push_back(glb_index);
          }
          else
          {
            cl.push_back(distdata[0].get_local(glb_index));
          }
        }

        // Mark cell for refinement
        if (cell_marker) cm.push_back(old_cell_marker->get(*c));
      }
    }

    // Add new cells
    uint cell_n = 0;
    for (uint i = 0; i < num_cells; i++)
    {
      if (cell_marker && cell_n == ndims)
      {
        cm.push_back(recv_buff_cell[i] == 1);
        cell_n = 0;
      }
      else
      {
        if (distdata[0].has_global(recv_buff_cell[i]))
        {
          cl.push_back(distdata[0].get_local(recv_buff_cell[i]));
        }
        else
        {
          cl.push_back(vi);
          for (uint j = 0; j < gdim; ++j)
          {
            coords.push_back(0.0);
          }
          distdata[0].set_map(vi, recv_buff_cell[i]);
          //FIXME: fix garbage
          //distdata.set_ghost(vi++, 0);
          shared_buffer.push_back(recv_buff_cell[i]);
        }
        cell_n++;
      }
    }

    // Exchange ghosted entities
    Array<real> send_buff;
    Array<uint> send_buff_indices, recv_source;
    send_size = shared_buffer.size();
    recv_count_vertices = static_cast<uint>(gdim) * send_size;
    recv_buff = new double[recv_count_vertices];
    double *rp = &recv_buff[0];
    recv_buff_map = new uint[send_size];
    uint *rmp = &recv_buff_map[0];
    MPI_Allreduce(&send_size, &recv_count, 1, MPI_INT, MPI_MAX,
                  MPI::DOLFIN_COMM);
    uint *shared = new uint[recv_count];
    for (uint i = 1; i < pe_size; i++)
    {

      src = (rank - i + pe_size) % pe_size;
      dest = (rank + i) % pe_size;

      MPI_Sendrecv(&shared_buffer[0], shared_buffer.size(), MPI_UNSIGNED, dest,
                   1, shared, recv_count, MPI_UNSIGNED, src, 1,
                   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recv_size);

      for (int j = 0; j < recv_size; j++)
        if (distdata[0].has_global(shared[j]))
        {
          if (!distdata[0].is_ghost(distdata[0].get_local(shared[j])))
          {
            offset = distdata[0].get_local(shared[j]) * gdim;
            for (uint d = 0; d < gdim; ++d)
            {
              send_buff.push_back(coords[offset + d]);
            }
            send_buff_indices.push_back(shared[j]);
            distdata[0].set_shared(distdata[0].get_local(shared[j]));
          }
          distdata[0].set_shared_adj(distdata[0].get_local(shared[j]), src);
        }

      MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_DOUBLE, src, 2, rp,
                   recv_count_vertices, MPI_DOUBLE, dest, 2, MPI::DOLFIN_COMM,
                   &status);
      MPI_Get_count(&status, MPI_DOUBLE, &recv_size);

      rp += recv_size;
      recv_count_vertices -= recv_size;

      MPI_Sendrecv(&send_buff_indices[0], send_buff_indices.size(),
                   MPI_UNSIGNED, src, 3, rmp, send_size, MPI_UNSIGNED, dest, 3,
                   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recv_size);

      rmp += recv_size;
      send_size -= recv_size;

      for (int k = 0; k < recv_size; k++)
        recv_source.push_back(status.MPI_SOURCE);

      send_buff.clear();
      send_buff_indices.clear();
    }

    uint j = 0;
    for (uint i = 0; i < shared_buffer.size(); i++)
    {
      offset = distdata[0].get_local(recv_buff_map[i]) * gdim;
      for (uint d = 0; d < gdim; ++d)
      {
        coords[offset + d] = recv_buff[j + d];
      }
      j += gdim;
      distdata[0].set_ghost(distdata[0].get_local(recv_buff_map[i]),
                               recv_source[i]);
    }
    shared_buffer.clear();
    recv_source.clear();
    delete[] recv_buff_map;
    delete[] recv_buff;
    delete[] shared;

  }
  delete[] recv_buff_cell;

  num_vertices = coords.size() / gdim;
  num_cells = cl.size() / ndims;

  // Construct new mesh and add all buffered entities
  Mesh new_mesh;
  MeshEditor editor(new_mesh, mesh.type(), mesh.geometry().dim());

//  for (uint d = 0; d < mesh.topology().dim(); ++d)
//  {
//    distdata.set_num_global(d, mesh.distdata().global_size(d));
//  }

  editor.init_vertices(num_vertices);
  editor.init_cells(num_cells);

  vi = 0;
  for (uint i = 0; i < coords.size(); i += gdim)
  {
    editor.add_vertex(vi++, &coords[i]);
  }
  coords.clear();

  dolfin_assert(cl.size() % ndims == 0);
  uint ci = 0;
  for (uint i = 0; i < cl.size(); i += ndims)
  {
    editor.add_cell(ci++, &cl[i]);
  }
  cl.clear();
  editor.close();

  // Overwrite old mesh with new, and invalidate numbering
  new_mesh.topology().distdata() = distdata;
  mesh = new_mesh;
  dolfin_assert(mesh.is_distributed());
//  mesh.distdata().set_invalid_numbering();
//  mesh.distdata().set_invalid_ownership();

  // Mark new cells for refinement
  if (cell_marker)
  {
    dolfin_assert(cm.size() == mesh.numCells());
    cell_marker->init(mesh, mesh.topology().dim());
    for (uint i = 0; i < cm.size(); i++)
    {
      cell_marker->set(i, cm[i]);
    }
    cm.clear();
  }
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distributeCommon(
    Mesh& mesh,
    MeshFunction<uint>& distribution,
    Array<std::pair<MeshFunction<uint> *, MeshFunction<uint> *> > *cell_functions,
    Array<std::pair<MeshFunction<real> *, MeshFunction<real> *> > *vertex_functions)
{
  typedef Array<std::pair<MeshFunction<uint> *, MeshFunction<uint> *> > CellFunctionArrayType;
  typedef Array<std::pair<MeshFunction<real> *, MeshFunction<real> *> > VertexFunctionArrayType;

  // new local distdata
  MeshDistributedData distdata(mesh.topology().dim());

  // sizes and numbers
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint gdim = mesh.geometry().dim();
  uint ndims = mesh.type().num_vertices(mesh.topology().dim());
  uint ncfunctions = (cell_functions) ? cell_functions->size() : 0;
  uint nvfunctions = (vertex_functions) ? vertex_functions->size() : 0;

  // vertex coordinates that have to be sent
  Array<real> *send_list_vertices = new Array<real> [pe_size];

  // global indices of vertices that have to be sent
  Array<uint> *send_list_mappings = new Array<uint> [pe_size];

  // indices of cell vertices (+ cell_function data) - to be sent
  Array<uint> *send_list_cells = new Array<uint> [pe_size];

  // vertex coordinates that are kept locally
  Array<real> coords;

  // cell list
  Array<uint> cl;

  // global indices of ghosted vertices
  Array<uint> shared_buffer;

  // exchanged function values
  Array<uint> *cfunctions = new Array<uint> [ncfunctions];
  Array<double> *vfunctions = new Array<double> [nvfunctions];

  uint num_cells, num_vertices, target_proc, glb_index, offset;

  // recv sizes and counters
  int recv_size, recv_size_cell, send_size;
  int recv_count, recv_count_vertices, recv_count_cells;

  //
  // Process mesh entities according to distribution
  //

  // counter for new local vertex indices
  uint vi = 0;

  // Distribution defined per vertex
  if (distribution.dim() == 0)
  {
    for (VertexIterator v(mesh); !v.end(); ++v)
    {
      if (!v->is_ghost())
      {
        // vertex has to be sent
        uint const glb_index = v->global_index();
        if (distribution.get(*v) != rank)
        {
          target_proc = distribution.get(*v);
          send_list_mappings[target_proc].push_back(glb_index);
          for (uint d = 0; d < gdim; ++d)
          {
            send_list_vertices[target_proc].push_back(v->x(d));
          }

          // Transfer vertex functions
          if (vertex_functions)
          {
            for (VertexFunctionArrayType::iterator f_it(
                vertex_functions->begin()); f_it != vertex_functions->end();
                ++f_it)
            {
              send_list_vertices[target_proc].push_back(f_it->first->get(*v));
            }
          }
        }
        // vertex stays here
        else if (distribution.get(*v) == rank)
        {
          for (uint d = 0; d < gdim; ++d)
          {
            coords.push_back(v->x(d));
          }
          distdata[0].set_map(vi++, glb_index);

          // store function values
          if (vertex_functions)
          {
            uint f_id(0);
            for (VertexFunctionArrayType::iterator f_it(
                vertex_functions->begin()); f_it != vertex_functions->end();
                ++f_it, ++f_id)
            {
              vfunctions[f_id].push_back(f_it->first->get(*v));
            }
          }
        }
      }
    }
    recv_count_cells = 0;
  }
  // Distribution defined per cell
  else if (distribution.dim() == mesh.topology().dim())
  {
    // mark vertices that have been added to send list already
    MeshFunction<bool> vertex_used(mesh, 0);
    vertex_used = false;

    for (CellIterator c(mesh); !c.end(); ++c)
    {
      // cell has to be sent
      if (distribution.get(*c) != rank)
      {
        target_proc = distribution.get(*c);
        for (VertexIterator v(*c); !v.end(); ++v)
        {
          // Buffer cell global vertex indices
          uint const glb_index = v->global_index();
          send_list_cells[target_proc].push_back(glb_index);

          // Buffer all cell vertices that belong to another processor
          if (!v->is_ghost() && !vertex_used.get(*v))
          {
            send_list_mappings[target_proc].push_back(glb_index);
            for (uint d = 0; d < gdim; ++d)
            {
              send_list_vertices[target_proc].push_back(v->x(d));
            }
            vertex_used.set(*v, true);

            // Transfer vertex functions
            if (vertex_functions)
            {
              for (VertexFunctionArrayType::iterator f_it(
                  vertex_functions->begin()); f_it != vertex_functions->end();
                  ++f_it)
              {
                send_list_vertices[target_proc].push_back(f_it->first->get(*v));
              }
            }
          }
        }

        // Transfer cell functions
        if (cell_functions)
        {
          for (CellFunctionArrayType::iterator f_it(cell_functions->begin());
              f_it != cell_functions->end(); ++f_it)
          {
            send_list_cells[target_proc].push_back(f_it->first->get(*c));
          }
        }
      }
      // cell stays here
      else
      {
        for (VertexIterator v(*c); !v.end(); ++v)
        {
          if (!vertex_used.get(*v))
          {
            if (!v->is_ghost())
            {
              for (uint d = 0; d < gdim; ++d)
              {
                coords.push_back(v->x(d));
              }
              distdata[0].set_map(vi++, v->global_index());
              vertex_used.set(*v, true);

              // store function values
              if (vertex_functions)
              {
                uint f_id(0);
                for (VertexFunctionArrayType::iterator f_it(
                    vertex_functions->begin()); f_it != vertex_functions->end();
                    ++f_it, ++f_id)
                {
                  vfunctions[f_id].push_back(f_it->first->get(*v));
                }
              }
            }
          }
        }
      }
    }

    // determine send/recv sizes for cells
    recv_count_cells = 0;
    for (uint i = 0; i < pe_size; i++)
    {
      send_size = send_list_cells[i].size();
      MPI_Reduce(&send_size, &recv_count_cells, 1, MPI_INT, MPI_SUM, i,
                 MPI::DOLFIN_COMM);
    }
  }
  // distribution defined on unknown entity
  else error("Distribution defined on unkown mesh entity");

  //
  // Exchange the processed entities
  //

  // determine send/recv sizes for vertices
  recv_count = 0;
  for (uint i = 0; i < pe_size; i++)
  {
    send_size = send_list_vertices[i].size();
    MPI_Reduce(&send_size, &recv_count, 1, MPI_INT, MPI_MAX, i,
               MPI::DOLFIN_COMM);
  }

  // number of vertices
  recv_count_vertices = recv_count / gdim;
  //num_vertices = recv_count;
  num_vertices = 0;

  // buffer for vertex coordinates
  double *recv_buff = new double[recv_count];

  // buffer for global vertex indices
  uint *recv_buff_map = new uint[recv_count_vertices];

  // number of cells
  num_cells = recv_count_cells;

  // buffer for cell data
  uint *recv_buff_cell = new uint[recv_count_cells];

  // set write pointer to beginning
  uint *rcp = &recv_buff_cell[0];

  // pairwise communication
  MPI_Status status;
  uint src, dest, buff_map;
  for (uint i = 1; i < pe_size; i++)
  {
    src = (rank - i + pe_size) % pe_size;
    dest = (rank + i) % pe_size;

    // exchange vertex coordinates
    MPI_Sendrecv(&send_list_vertices[dest][0], send_list_vertices[dest].size(),
                 MPI_DOUBLE, dest, 0, recv_buff, recv_count, MPI_DOUBLE, src, 0,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &recv_size);

    // exchange global indices
    MPI_Sendrecv(&send_list_mappings[dest][0], send_list_mappings[dest].size(),
                 MPI_UNSIGNED, dest, 1, recv_buff_map, recv_count_vertices,
                 MPI_UNSIGNED, src, 1, MPI::DOLFIN_COMM, &status);

    // exchange cells
    MPI_Sendrecv(&send_list_cells[dest][0], send_list_cells[dest].size(),
                 MPI_UNSIGNED, dest, 2, rcp, recv_count_cells, MPI_UNSIGNED,
                 src, 2, MPI::DOLFIN_COMM, &status);

    // move write pointer
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_size_cell);
    rcp += recv_size_cell;
    recv_count_cells -= recv_size_cell;

    // process vertices and build map of local <-> global indices
    buff_map = 0;
    for (int i = 0; i < recv_size; i += gdim + nvfunctions)
    {
      if (!distdata[0].has_global(recv_buff_map[buff_map]))
      {
        distdata[0].set_map(vi++, recv_buff_map[buff_map]);
        for (uint d = 0; d < gdim; ++d)
        {
          coords.push_back(recv_buff[i + d]);
        }

        // store function values
        if (vertex_functions)
        {
          for (uint f_id(0); f_id < nvfunctions; ++f_id)
          {
            vfunctions[f_id].push_back(recv_buff[i + gdim + f_id]);
          }
        }
      }
      buff_map++;
    }
  }

  //
  // Clear send buffers
  //

  for (uint i = 0; i < pe_size; i++)
  {
    send_list_cells[i].clear();
    send_list_vertices[i].clear();
    send_list_mappings[i].clear();
  }
  delete[] send_list_vertices;
  delete[] send_list_mappings;
  delete[] send_list_cells;
  delete[] recv_buff_map;
  delete[] recv_buff;

  //
  // Process cells if distribution is defined on cells
  //

  if (distribution.dim() == mesh.topology().dim())
  {
    // Add cells that stayed
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      if (distribution.get(*c) == rank)
      {
        for (VertexIterator v(*c); !v.end(); ++v)
        {
          uint const glb_index = v->global_index();
          // ghosted vertex
          if (!distdata[0].has_global(glb_index))
          {
            cl.push_back(vi);
            for (uint j = 0; j < gdim; j++)
              coords.push_back(0.0);
            distdata[0].set_map(vi, glb_index);
            //distdata.set_ghost(vi++, 0);
            shared_buffer.push_back(glb_index);

            if (vertex_functions)
            {
              for (uint f_id(0); f_id < nvfunctions; ++f_id)
              {
                vfunctions[f_id].push_back(0);
              }
            }
          }
          // own vertex
          else
          {
            cl.push_back(distdata[0].get_local(glb_index));
          }
        }

        // Conserve function values
        if (cell_functions)
        {
          uint f_id(0);
          for (CellFunctionArrayType::iterator f_it(cell_functions->begin());
              f_it != cell_functions->end(); ++f_it, ++f_id)
          {
            cfunctions[f_id].push_back(f_it->first->get(*c));
          }
        }
      }
    }

    // Add cells that have been received
    uint cell_n = 0;
    for (uint i = 0; i < num_cells; i++)
    {
      // function values
      if (cell_functions && cell_n >= ndims)
      {
        cfunctions[cell_n - ndims].push_back(recv_buff_cell[i]);
        ++cell_n;
        if ((cell_n - ndims) == ncfunctions) cell_n = 0;
      }
      // vertex indices
      else
      {
        // own vertex
        if (distdata[0].has_global(recv_buff_cell[i])) cl.push_back(
            distdata[0].get_local(recv_buff_cell[i]));
        // ghosted vertex
        else
        {
          cl.push_back(vi);
          for (uint j = 0; j < gdim; ++j)
          {
            coords.push_back(0.0);
          }
          distdata[0].set_map(vi, recv_buff_cell[i]);
          //distdata.set_ghost(vi++, 0);
          shared_buffer.push_back(recv_buff_cell[i]);

          if (vertex_functions)
          {
            for (uint f_id(0); f_id < nvfunctions; ++f_id)
            {
              //vfunctions[f_id].push_back(recv_buff[i+gdim+f_id]);
              vfunctions[f_id].push_back(0);
            }
          }
        }
        cell_n++;
      }
    }

    //
    // Exchange ghosted entities
    //

    // ghost vertex coordinates, that are going to be sent
    Array<real> send_buff;

    // ghosted vertex indices
    Array<uint> send_buff_indices, recv_source;

    // number of requested ghost vertices
    send_size = shared_buffer.size();
    recv_count_vertices = (static_cast<uint>(gdim) + nvfunctions) * send_size;

    // buffer for received ghost coordinates
    recv_buff = new double[recv_count_vertices];
    double *rp = &recv_buff[0];

    // buffer for received global ids of ghosts
    recv_buff_map = new uint[send_size];
    uint *rmp = &recv_buff_map[0];

    // determine maximum message size
    MPI_Allreduce(&send_size, &recv_count, 1, MPI_INT, MPI_MAX,
                  MPI::DOLFIN_COMM);

    // buffer for requested indices
    uint *shared = new uint[recv_count];

    // pairwise communication
    for (uint i = 1; i < pe_size; i++)
    {

      src = (rank - i + pe_size) % pe_size;
      dest = (rank + i) % pe_size;

      // exchange indices for requested vertices
      MPI_Sendrecv(&shared_buffer[0], shared_buffer.size(), MPI_UNSIGNED, dest,
                   1, shared, recv_count, MPI_UNSIGNED, src, 1,
                   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recv_size);

      // prepare coordinates and global indices of requested vertices
      for (int j = 0; j < recv_size; j++)
      {
        if (distdata[0].has_global(shared[j]))
        {
          if (!distdata[0].is_ghost(distdata[0].get_local(shared[j])))
          {
            offset = distdata[0].get_local(shared[j]) * gdim;
            for (uint d = 0; d < gdim; ++d)
            {
              send_buff.push_back(coords[offset + d]);
            }

            if (vertex_functions)
            {
              for (uint f_id(0); f_id < nvfunctions; ++f_id)
              {
                send_buff.push_back(
                    vfunctions[f_id][distdata[0].get_local(shared[j])]);
              }
            }

            send_buff_indices.push_back(shared[j]);
            distdata[0].set_shared(distdata[0].get_local(shared[j]));
          }
          distdata[0].set_shared_adj(distdata[0].get_local(shared[j]), src);
        }
      }

      // exchange prepared coordinates
      MPI_Sendrecv(&send_buff[0], send_buff.size(), MPI_DOUBLE, src, 2, rp,
                   recv_count_vertices, MPI_DOUBLE, dest, 2, MPI::DOLFIN_COMM,
                   &status);
      MPI_Get_count(&status, MPI_DOUBLE, &recv_size);

      // move write pointer
      rp += recv_size;
      recv_count_vertices -= recv_size;

      // exchange prepared indices
      MPI_Sendrecv(&send_buff_indices[0], send_buff_indices.size(),
                   MPI_UNSIGNED, src, 3, rmp, send_size, MPI_UNSIGNED, dest, 3,
                   MPI::DOLFIN_COMM, &status);
      MPI_Get_count(&status, MPI_UNSIGNED, &recv_size);

      // move write pointer
      rmp += recv_size;
      send_size -= recv_size;

      for (int k = 0; k < recv_size; k++)
      {
        recv_source.push_back(status.MPI_SOURCE);
      }

      send_buff.clear();
      send_buff_indices.clear();
    }

    // update coordinates of received ghost vertices
    uint j = 0;
    for (uint i = 0; i < shared_buffer.size(); i++)
    {
      offset = distdata[0].get_local(recv_buff_map[i]) * gdim;
      for (uint d = 0; d < gdim; ++d)
      {
        coords[offset + d] = recv_buff[j + d];
      }
      j += gdim;
      distdata[0].set_ghost(distdata[0].get_local(recv_buff_map[i]),
                               recv_source[i]);
      if (vertex_functions)
      {
        for (uint f_id(0); f_id < nvfunctions; ++f_id)
        {
          vfunctions[f_id][distdata[0].get_local(recv_buff_map[i])] =
              recv_buff[j + f_id];
        }
        j += nvfunctions;
      }
    }
    shared_buffer.clear();
    recv_source.clear();
    delete[] recv_buff_map;
    delete[] recv_buff;
    delete[] shared;

  }
  delete[] recv_buff_cell;

  //
  // Build the new local mesh
  //

  // new cell and vertex numbers
  num_vertices = coords.size() / gdim;
  num_cells = cl.size() / ndims;

  // Construct new mesh and add all buffered entities
  Mesh new_mesh;
  MeshEditor editor(new_mesh, mesh.type(), mesh.geometry().dim());

//  for (uint d = 0; d < mesh.topology().dim(); ++d)
//  {
//    distdata.set_num_global(d, mesh.distdata().global_size(d));
//  }

  editor.init_vertices(num_vertices);
  editor.init_cells(num_cells);

  // Add vertices
  vi = 0;
  for (uint i = 0; i < coords.size(); i += gdim)
  {
    editor.add_vertex(vi++, &coords[i]);
  }
  coords.clear();

  // Add cells
  dolfin_assert(cl.size() % ndims == 0);
  uint ci = 0;
  for (uint i = 0; i < cl.size(); i += ndims)
  {
    editor.add_cell(ci++, &cl[i]);
  }
  cl.clear();
  editor.close();

  // Overwrite old mesh with new, and invalidate numbering
  new_mesh.topology().distdata() = distdata;
  mesh = new_mesh;
  dolfin_assert(mesh.is_distributed());
//  mesh.distdata().set_invalid_numbering();
//  mesh.distdata().set_invalid_ownership();

  // Fill new cell functions
  if (cell_functions && mesh.numCells() > 0)
  {
    uint f_id(0);
    for (CellFunctionArrayType::iterator f_it(cell_functions->begin());
        f_it != cell_functions->end(); ++f_it, ++f_id)
    {
      dolfin_assert(cfunctions[f_id].size() == mesh.numCells());
      f_it->second->init(mesh, mesh.topology().dim());
      for (uint i(0); i < cfunctions[f_id].size(); ++i)
      {
        f_it->second->set(i, cfunctions[f_id][i]);
      }
      cfunctions[f_id].clear();
    }
  }

  // Fill new vertex functions
  if (vertex_functions && mesh.numVertices() > 0)
  {
    uint f_id(0);
    for (VertexFunctionArrayType::iterator f_it(vertex_functions->begin());
        f_it != vertex_functions->end(); ++f_it, ++f_id)
    {
      dolfin_assert(vfunctions[f_id].size() == mesh.numVertices());
      f_it->second->init(mesh, 0);
      for (uint i(0); i < vfunctions[f_id].size(); ++i)
      {
        f_it->second->set(i, vfunctions[f_id][i]);
      }
      vfunctions[f_id].clear();
    }
  }

  // Delete function arrays
  delete[] cfunctions;
  delete[] vfunctions;
}
//-----------------------------------------------------------------------------

#else

//-----------------------------------------------------------------------------
MPIMeshCommunicator::MPIMeshCommunicator()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MPIMeshCommunicator::~MPIMeshCommunicator()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(Mesh& mesh,
    MeshFunction<uint>& distribution)
{

  error("Cannot distribute mesh without MPI.");
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(Mesh& mesh,
    MeshFunction<uint>& distribution,
    MeshFunction<bool>& old_cell_marker,
    MeshFunction<bool>& cell_marker)
{
  error("Cannot distribute mesh without MPI.");
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(Mesh& mesh,
    MeshFunction<uint>& distribution,
    Array< std::pair< MeshFunction<uint> *,
    MeshFunction<uint> * > >& cell_functions)
{
  error("Cannot distribute mesh without MPI.");
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(
    Mesh& mesh,
    MeshFunction<uint>& distribution,
    Array<std::pair<MeshFunction<real> *, MeshFunction<real> *> >& vertex_functions)
{
}
//-----------------------------------------------------------------------------
void MPIMeshCommunicator::distribute(
    Mesh& mesh,
    MeshFunction<uint>& distribution,
    Array<std::pair<MeshFunction<uint> *, MeshFunction<uint> *> >& cell_functions,
    Array<std::pair<MeshFunction<real> *, MeshFunction<real> *> >& vertex_functions)
{
  error("Cannot distribute mesh without MPI.");
}
//-----------------------------------------------------------------------------

#endif
