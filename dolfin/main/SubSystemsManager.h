// Copyright (C) 2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson 2009.
//
// First added:  2008-01-07
// Last changed: 2009-04-22

#ifndef __SUB_SYSTEMS_MANAGER_H
#define __SUB_SYSTEMS_MANAGER_H

namespace dolfin
{
  
  /// This is a singleton class which manages the initialisation and 
  /// finalisation of various sub systems, such as MPI and PETSc.

  class SubSystemsManager
  {
  public:

    /// Initialise MPI
    static void initMPI();

    /// Initialize PETSc without command-line arguments
    static void initPETSc();

    /// Initialize PETSc with command-line arguments
    static void initPETSc(int argc, char* argv[], bool cmd_line_args = true);

    /// Finalize MPI
    static void finalizeMPI();

    /// Finalize PETSc
    static void finalizePETSc();

  private:

    // Constructor
    SubSystemsManager();

    // Copy construtor
    SubSystemsManager(const SubSystemsManager& sub_sys_manager);

    // Destructor
    ~SubSystemsManager();

    // Check if MPI has been initialised (returns true if MPI has been 
    //   initialised, even if it is later finalised) 
    static bool MPIinitialized();

    // Singleton instance
    static SubSystemsManager sub_systems_manager;

    // State variables
    bool petsc_initialized;
    bool petsc_controls_mpi;

  };

}

#endif
