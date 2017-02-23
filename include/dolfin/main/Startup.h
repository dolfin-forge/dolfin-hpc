#ifndef __DOLFIN_STARTUP_H
#define __DOLFIN_STARTUP_H

#include <dolfin/main/init.h>

namespace dolfin
{

struct Startup
{

  Startup(int argc, char* argv[]) { dolfin_init(argc, argv); }

  virtual ~Startup() { dolfin_finalize(); }

};

} /* namespace dolfin */

#endif /* __DOLFIN_STARTUP_H */
