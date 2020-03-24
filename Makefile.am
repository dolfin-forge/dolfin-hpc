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

EXTRA_DIST = data/meshes/aneurysm.bim\
             data/meshes/README\
             data/meshes/snake_hole.bim\
             data/meshes/circle.bim\
             data/meshes/cylinder.bim\
             data/meshes/dolfin-1.bim\
             data/meshes/dolfin-2.bim\
             data/meshes/gear.bim\
             data/meshes/meshbc.bim\
             data/meshes/slit.bim\
             data/meshes/snake.bim\
             data/meshes/tetmesh-1.bim\
             data/meshes/sphere.bim\
             data/meshes/tetmesh_24_8_8.bim\
             data/meshes/tetmesh-16.bim\
             data/meshes/tetmesh-4.bim\
             data/meshes/tetmesh-8.bim\
             data/meshes/tetmesh_backward_facing_step_32_8_8.bim\
             data/meshes/tetmesh_LES_10_2_7.bim\
             data/meshes/tetmesh_bluff_body_32_8_8.bim\
             data/meshes/tetrahedron.bim\
             data/meshes/torso.bim\
             data/meshes/triangle.bim\
             data/meshes/trimesh-1.bim\
             data/meshes/trimesh-128.bim\
             data/meshes/trimesh-16.bim\
             data/meshes/trimesh-2.bim\
             data/meshes/trimesh-32.bim\
             data/meshes/trimesh-4.bim\
             data/meshes/trimesh-64.bim\
             data/meshes/trimesh-8.bim\
             data/meshes/two_particles.bim\
             data/meshes/unstruct-trimesh-1.bim\
             data/meshes/unstruct-trimesh-2.bim\
             demo/function/Makefile\
             demo/function/Projection.ufl\
             demo/function/Projection.h\
             demo/function/main.cpp\
             demo/insitu/poisson/Makefile\
             demo/insitu/poisson/main.cpp\
             demo/insitu/poisson/Poisson.ufl\
             demo/insitu/poisson/Poisson.h\
             demo/insitu/poisson/UnitSquareMesh_32x32.bim\
             demo/insitu/poisson/poisson.session\
             demo/insitu/convection-diffusion/Makefile\
             demo/insitu/convection-diffusion/main.cpp\
             demo/insitu/convection-diffusion/ConvectionDiffusion.ufl\
             demo/insitu/convection-diffusion/ConvectionDiffusion.h\
             demo/insitu/convection-diffusion/UnitSquareMesh_32x32.bim\
             demo/mesh/refinement/Makefile\
             demo/mesh/refinement/main.cpp\
             demo/mesh/meshfunction/mesh2D.bim\
             demo/mesh/meshfunction/meshfunction.bim\
             demo/mesh/meshfunction/Makefile\
             demo/mesh/meshfunction/main.cpp\
             demo/mesh/subdomains/Makefile\
             demo/mesh/subdomains/main.cpp\
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
             demo/pde/poisson/README\
             demo/pde/poisson/Makefile\
             demo/pde/poisson/main.cpp\
             demo/pde/poisson/Poisson.ufl\
             demo/pde/poisson/Poisson.h\
             demo/pde/stokes/taylor-hood/Makefile\
             demo/pde/stokes/taylor-hood/Stokes.ufl\
             demo/pde/stokes/taylor-hood/Stokes.h\
             demo/pde/stokes/taylor-hood/main.cpp\
             demo/quadrature/main.cpp\
             demo/quadrature/Makefile
