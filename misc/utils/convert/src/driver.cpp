/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <njansson@kth.se> wrote this file. As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return Niclas Jansson
 * ----------------------------------------------------------------------------
 */

#include <cstring>
#include <iostream>
#include <unistd.h>

#include <Mesh.h>
#include <DOLFINxml.h>
#include <DOLFINBinary.h>

int main(int argc, char * argv [])
{

  std::cout << "This is " PACKAGE_NAME <<
    " Version " << PACKAGE_VERSION << std::endl;
  if (argc < 3) {
    std::cerr << "\nUsage: dolfin-convert <input-mesh> <output-mesh>\n";
    return 1;
  }
  
  std::string input_fname, output_fname;
  input_fname.append(argv[1]);
  output_fname.append(argv[2]);

  Mesh *mesh;
  
  uint32_t const len = input_fname.size();
  
  // Choose file type base on suffix
  if (((input_fname.rfind(".xml") != std::string::npos) &&
       (input_fname.rfind(".xml") + std::string(".xml").size() == len)) ||
      ((input_fname.rfind(".xml.gz") != std::string::npos) &&
       (input_fname.rfind(".xml.gz") + std::string(".xml.gz").size() == len))) {
    mesh = new DOLFINxml();
  }
  else
  {
    std::cerr << "Unkown mesh type\n";
    return -1;
  }

  mesh->load_mesh(input_fname);

  write_dolfin_binary(mesh, output_fname);
  
  return 0;
}


