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

#include <dolfin/mesh/MPIMeshCommunicator.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

namespace dolfin
{

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
  if (!mesh.is_distributed())
  {
    return;
  }

#if HAVE_MPI

  message("MPIMeshCommunicator::distribute");
  mesh.distdata()[0].disp();

  MeshDistributedData distdata(mesh.topology().dim());
  uint rank = MPI::processNumber();
  uint pe_size = MPI::numProcesses();
  uint gdim = mesh.geometry().dim();
  uint ndims = mesh.type().num_entities(0);

  Array<real> *send_list_vertices = new Array<real> [pe_size];
  Array<uint> *send_list_mappings = new Array<uint> [pe_size];
  Array<uint> *send_list_cells = new Array<uint> [pe_size];
  Array<real> coords;
  Array<uint> cl, shared_buffer;
  Array<bool> cm;
  uint num_cells, num_vertices, target_proc, offset;

  int recv_size, recv_size_cell, send_size;
  int recv_count, recv_count_vertices, recv_count_cells;

  message("Process mesh entities");
  // Process mesh entities according to distribution
  uint vi = 0;
  // Distribution defined per vertex
  if (distribution.dim() == 0)
  {
    message(1, "MPIMeshCommunicator : vertex distribution");
    for (VertexIterator v(mesh); !v.end(); ++v)
    {
      if (v->is_owned())
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
    message(1, "MPIMeshCommunicator : cell distribution");
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
    error("Distribution defined on unknown mesh entity");
  }

  message("Exchange the processed entities");
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
    MPI_DOUBLE,
                 dest, 0, recv_buff, recv_count, MPI_DOUBLE, src, 0,
                 MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &recv_size);

    MPI_Sendrecv(&send_list_mappings[dest][0], send_list_mappings[dest].size(),
    MPI_UNSIGNED,
                 dest, 1, recv_buff_map, recv_count_vertices,
                 MPI_UNSIGNED,
                 src, 1, MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&send_list_cells[dest][0], send_list_cells[dest].size(),
    MPI_UNSIGNED,
                 dest, 2, rcp, recv_count_cells, MPI_UNSIGNED, src, 2,
                 MPI::DOLFIN_COMM, &status);
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
    message("Process new and old cells");
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
            message("%8u owned by %8u",vi, distribution.get(*c));
            distdata[0].set_ghost(vi++, distribution.get(*c));
            shared_buffer.push_back(glb_index);
          }
          else
          {
            cl.push_back(distdata[0].get_local(glb_index));
          }
        }
      }
    }

    // Add new cells
    uint cell_n = 0;
    for (uint i = 0; i < num_cells; i++)
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
        distdata[0].set_ghost(vi++, pe_size);
        shared_buffer.push_back(recv_buff_cell[i]);
      }
      cell_n++;
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
      MPI_UNSIGNED,
                   src, 3, rmp, send_size, MPI_UNSIGNED, dest, 3,
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
  message("Construct new mesh");
  Mesh new_mesh;
  MeshEditor editor(new_mesh, mesh.type(), mesh.geometry().dim());

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

  message("MPIMeshCommunicator::done");

#endif
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
