ACLOCAL_AMFLAGS = -I m4 --install
SUBDIRS = src include

if ENABLE_CHECK 
SUBDIRS += tests
endif

if ENABLE_PYTHON_MODULES
SUBDIRS += site-packages
endif

if ENABLE_DOXYGEN
SUBDIRS += doc/progr_reference
endif

dist_data_DATA = regen.sh

EXTRA_DIST = data/meshes/aneurysm.xml.gz\
             data/meshes/README\
             data/meshes/snake_hole.xml.gz\
             data/meshes/circle.xml.gz\
             data/meshes/cylinder.xml.gz\
             data/meshes/dolfin-1.xml.gz\
             data/meshes/dolfin-2.xml.gz\
             data/meshes/gear.xml.gz\
             data/meshes/meshbc.xml.gz\
             data/meshes/slit.xml.gz\
             data/meshes/snake.xml.gz\
             data/meshes/tetmesh-1.xml.gz\
             data/meshes/sphere.xml.gz\
             data/meshes/tetmesh_24_8_8.xml.gz\
             data/meshes/tetmesh-16.xml.gz\
             data/meshes/tetmesh-4.xml.gz\
             data/meshes/tetmesh-8.xml.gz\
             data/meshes/tetmesh_backward_facing_step_32_8_8.xml.gz\
             data/meshes/tetmesh_LES_10_2_7.xml.gz\
             data/meshes/tetmesh_bluff_body_32_8_8.xml.gz\
             data/meshes/tetrahedron.xml.gz\
             data/meshes/torso.xml.gz\
             data/meshes/triangle.xml.gz\
             data/meshes/trimesh-1.xml.gz\
             data/meshes/trimesh-128.xml.gz\
             data/meshes/trimesh-16.xml.gz\
             data/meshes/trimesh-2.xml.gz\
             data/meshes/trimesh-32.xml.gz\
             data/meshes/trimesh-4.xml.gz\
             data/meshes/trimesh-64.xml.gz\
             data/meshes/trimesh-8.xml.gz\
             data/meshes/two_particles.xml.gz\
             data/meshes/unstruct-trimesh-1.xml.gz\
             data/meshes/unstruct-trimesh-2.xml.gz\
             demo/function/Makefile\
             demo/function/Projection.ufl\
             demo/function/Projection.h\
             demo/function/main.cpp\
             demo/insitu/poisson/Makefile\
             demo/insitu/poisson/main.cpp\
             demo/insitu/poisson/Poisson.ufl\
             demo/insitu/poisson/Poisson.h\
	     demo/insitu/poisson/UnitSquareMesh_32x32.xml\
	     demo/insitu/poisson/poisson.session\
             demo/insitu/convection-diffusion/Makefile\
             demo/insitu/convection-diffusion/main.cpp\
             demo/insitu/convection-diffusion/ConvectionDiffusion.ufl\
             demo/insitu/convection-diffusion/ConvectionDiffusion.h\
	     demo/insitu/convection-diffusion/UnitSquareMesh_32x32.xml\
             demo/mesh/ale/Makefile\
             demo/mesh/ale/main.cpp\
             demo/mesh/refinement/Makefile\
             demo/mesh/refinement/main.cpp\
             demo/mesh/meshfunction/mesh2D.xml.gz\
             demo/mesh/meshfunction/meshfunction.xml\
             demo/mesh/meshfunction/Makefile\
             demo/mesh/meshfunction/main.cpp\
             demo/mesh/partitioning/main.cpp\
             demo/mesh/partitioning/Makefile\
             demo/mesh/subdomains/Makefile\
             demo/mesh/subdomains/main.cpp\
             demo/nls/cahn-hilliard/CahnHilliard2D.ufl\
             demo/nls/cahn-hilliard/CahnHilliard2D.h\
             demo/nls/cahn-hilliard/CahnHilliard3D.ufl\
             demo/nls/cahn-hilliard/CahnHilliard3D.h\
             demo/nls/cahn-hilliard/Makefile\
             demo/nls/cahn-hilliard/main.cpp\
             demo/nls/nonlinearpoisson/Makefile\
             demo/nls/nonlinearpoisson/NonlinearPoisson.ufl\
             demo/nls/nonlinearpoisson/NonlinearPoisson.h\
             demo/nls/nonlinearpoisson/main.cpp\
             demo/pde/convection-diffusion/Makefile\
             demo/pde/convection-diffusion/ConvectionDiffusion.ufl\
             demo/pde/convection-diffusion/subdomains.xml.gz\
             demo/pde/convection-diffusion/velocity.xml.gz\
             demo/pde/convection-diffusion/ConvectionDiffusion.h\
             demo/pde/convection-diffusion/main.cpp\
             demo/pde/dg/advection_diffusion/main.cpp\
             demo/pde/dg/advection_diffusion/AdvectionDiffusion.ufl\
             demo/pde/dg/advection_diffusion/velocity.xml.gz\
             demo/pde/dg/advection_diffusion/AdvectionDiffusion.h\
             demo/pde/dg/advection_diffusion/Makefile\
             demo/pde/dg/advection_diffusion/OutflowFacet.ufl\
             demo/pde/dg/advection_diffusion/OutflowFacet.h\
             demo/pde/dg/advection_diffusion/Projection.ufl\
             demo/pde/dg/advection_diffusion/Projection.h\
             demo/pde/dg/poisson/main.cpp\
             demo/pde/dg/poisson/Makefile\
             demo/pde/dg/poisson/P1Projection.ufl\
             demo/pde/dg/poisson/P1Projection.h\
             demo/pde/dg/poisson/Poisson.ufl\
             demo/pde/dg/poisson/Poisson.h\
             demo/pde/elasticity/main.cpp\
             demo/pde/elasticity/Elasticity.ufl\
             demo/pde/elasticity/Makefile\
             demo/pde/elasticity/Elasticity.h\
             demo/pde/functional/EnergyNorm.h\
             demo/pde/functional/EnergyNorm.ufl\
             demo/pde/functional/Makefile\
             demo/pde/functional/main.cpp\
             demo/pde/lift-drag/Drag.ufl\
             demo/pde/lift-drag/pressure.xml.gz\
             demo/pde/lift-drag/Drag.h\
             demo/pde/lift-drag/Lift.ufl\
             demo/pde/lift-drag/Lift.h\
             demo/pde/lift-drag/Makefile\
             demo/pde/lift-drag/main.cpp\
             demo/pde/mixed-poisson/MixedPoisson.ufl\
             demo/pde/mixed-poisson/MixedPoisson.h\
             demo/pde/mixed-poisson/P1Projection.ufl\
             demo/pde/mixed-poisson/P1Projection.h\
             demo/pde/mixed-poisson/Makefile\
             demo/pde/mixed-poisson/main.cpp\
             demo/pde/nonlinear-poisson/Makefile\
             demo/pde/nonlinear-poisson/NonlinearPoisson.ufl\
             demo/pde/nonlinear-poisson/NonlinearPoisson.h\
             demo/pde/nonlinear-poisson/main.cpp\
             demo/pde/poisson/README\
             demo/pde/poisson/Makefile\
             demo/pde/poisson/main.cpp\
             demo/pde/poisson/Poisson.ufl\
             demo/pde/poisson/Poisson.h\
             demo/pde/stokes/Stokes2D.h\
             demo/pde/stokes/stabilized/Makefile\
             demo/pde/stokes/stabilized/Stokes.ufl\
             demo/pde/stokes/stabilized/Stokes.h\
             demo/pde/stokes/stabilized/main.cpp\
             demo/pde/stokes/taylor-hood/Makefile\
             demo/pde/stokes/taylor-hood/Stokes.ufl\
             demo/pde/stokes/taylor-hood/Stokes.h\
             demo/pde/stokes/taylor-hood/main.cpp\
             demo/quadrature/main.cpp\
             demo/quadrature/Makefile
