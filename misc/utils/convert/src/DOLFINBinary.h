/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <njansson@kth.se> wrote this file. As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return Niclas Jansson
 * ----------------------------------------------------------------------------
 */

#ifndef DOLFIN_CONVERT_BINARY_H
#define DOLFIN_CONVERT_BINARY_H

#define BINARY_MAGIC 0xB4B3

#include <string>
#include <fstream>

enum Binary_data_t
{
  BINARY_MESH_DATA,
  BINARY_VECTOR_DATA,
  BINARY_FUNCTION_DATA,
  BINARY_MESH_FUNCTION_DATA
};

typedef struct
{
  uint32_t magic;
  uint32_t bendian;
  uint32_t pe_size;
  Binary_data_t type;
} BinaryFileHeader;


void write_dolfin_binary(Mesh *mesh, std::string& filename)
{
  std::ofstream fp(filename.c_str(), std::ofstream::binary);
  if ( not fp.good() )
    throw std::runtime_error( "Failed to open file \"" + filename + "\"" );

  BinaryFileHeader hdr;
  hdr.magic = BINARY_MAGIC;
  
  /* Write Header */
  fp.write(reinterpret_cast<const char*>(&hdr), sizeof(BinaryFileHeader));
  fp.write(reinterpret_cast<const char*>(&mesh->gdim), sizeof(uint32_t));
  fp.write(reinterpret_cast<const char*>(&mesh->cell_type), sizeof(uint32_t));

  /* Write vertices */
  fp.write(reinterpret_cast<const char*>(&mesh->num_vertices), sizeof(uint32_t));
  fp.write(reinterpret_cast<const char*>(mesh->vertices),
	   mesh->gdim * mesh->num_vertices * sizeof(double));

  /* Write cells */
  fp.write(reinterpret_cast<const char*>(&mesh->num_cells), sizeof(uint32_t));
  fp.write(reinterpret_cast<const char*>(mesh->cells),
	   mesh->num_cells * (mesh->cell_type + 1) * sizeof(int));

  fp.close();

  std::cout << "DOLFIN Binary written to \"" + filename + "\"\n";
  
}

#endif
