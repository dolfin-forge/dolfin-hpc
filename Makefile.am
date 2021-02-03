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

if ENABLE_MESH_CONVERTER
SUBDIRS += misc/utils/convert
endif

dist_data_DATA = regen.sh

EXTRA_DIST = data/meshes/aneurysm.bin\
             data/meshes/README\
             data/meshes/snake_hole.bin\
             data/meshes/circle.bin\
             data/meshes/cylinder.bin\
             data/meshes/dolfin-1.bin\
             data/meshes/dolfin-2.bin\
             data/meshes/gear.bin\
             data/meshes/meshbc.bin\
             data/meshes/slit.bin\
             data/meshes/snake.bin\
             data/meshes/tetmesh-1.bin\
             data/meshes/sphere.bin\
             data/meshes/tetmesh_24_8_8.bin\
             data/meshes/tetmesh-16.bin\
             data/meshes/tetmesh-4.bin\
             data/meshes/tetmesh-8.bin\
             data/meshes/tetmesh_backward_facing_step_32_8_8.bin\
             data/meshes/tetmesh_LES_10_2_7.bin\
             data/meshes/tetmesh_bluff_body_32_8_8.bin\
             data/meshes/tetrahedron.bin\
             data/meshes/torso.bin\
             data/meshes/triangle.bin\
             data/meshes/trimesh-1.bin\
             data/meshes/trimesh-128.bin\
             data/meshes/trimesh-16.bin\
             data/meshes/trimesh-2.bin\
             data/meshes/trimesh-32.bin\
             data/meshes/trimesh-4.bin\
             data/meshes/trimesh-64.bin\
             data/meshes/trimesh-8.bin\
             data/meshes/two_particles.bin\
             data/meshes/unstruct-trimesh-1.bin\
             data/meshes/unstruct-trimesh-2.bin\
             demo/function/Makefile\
             demo/function/Projection.ufl\
             demo/function/Projection.h\
             demo/function/main.cpp\
             demo/insitu/poisson/Makefile\
             demo/insitu/poisson/main.cpp\
             demo/insitu/poisson/Poisson.ufl\
             demo/insitu/poisson/Poisson.h\
             demo/insitu/poisson/UnitSquareMesh_32x32.bin\
             demo/insitu/poisson/poisson.session\
             demo/insitu/convection-diffusion/Makefile\
             demo/insitu/convection-diffusion/main.cpp\
             demo/insitu/convection-diffusion/ConvectionDiffusion.ufl\
             demo/insitu/convection-diffusion/ConvectionDiffusion.h\
             demo/insitu/convection-diffusion/UnitSquareMesh_32x32.bin\
             demo/mesh/refinement/Makefile\
             demo/mesh/refinement/main.cpp\
             demo/mesh/meshfunction/mesh2D.bin\
             demo/mesh/meshfunction/meshfunction.bin\
             demo/mesh/meshfunction/Makefile\
             demo/mesh/meshfunction/main.cpp\
             demo/mesh/subdomains/Makefile\
             demo/mesh/subdomains/main.cpp\
             demo/pde/bcs/main.cpp\
             demo/pde/bcs/Makefile\
             demo/pde/bcs/Poisson.ufl\
             demo/pde/bcs/Poisson.h\
             demo/pde/dg/main.cpp\
             demo/pde/dg/Makefile\
             demo/pde/dg/P1Projection.ufl\
             demo/pde/dg/P1Projection.h\
             demo/pde/dg/Poisson.ufl\
             demo/pde/dg/Poisson.h\
             demo/pde/elasticity/main.cpp\
             demo/pde/elasticity/Elasticity.ufl\
             demo/pde/elasticity/Makefile\
             demo/pde/elasticity/Elasticity.h\
             demo/pde/functional/EnergyNorm.h\
             demo/pde/functional/EnergyNorm.ufl\
             demo/pde/functional/Makefile\
             demo/pde/functional/main.cpp\
             demo/pde/poisson/README\
             demo/pde/poisson/Makefile\
             demo/pde/poisson/main.cpp\
             demo/pde/poisson/Poisson.ufl\
             demo/pde/poisson/Poisson.h\
             demo/pde/poisson/UnitSquareMesh_32x32.bin\
             demo/pde/stokes/Stokes2D.h\
             demo/pde/stokes/Makefile\
             demo/pde/stokes/Stokes.ufl\
             demo/pde/stokes/Stokes.h\
             demo/pde/stokes/main.cpp\
             demo/quadrature/main.cpp\
             demo/quadrature/Makefile
