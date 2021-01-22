SUBDIRS=celltypes entities mappings unitmeshes

AM_CPPFLAGS = -I@top_srcdir@/include
noinst_LTLIBRARIES = libdolfin_mesh.la

libdolfin_mesh_la_LIBADD = celltypes/libdolfin_mesh_celltypes.la\
                           entities/libdolfin_mesh_entities.la\
                           mappings/libdolfin_mesh_mappings.la\
                           unitmeshes/libdolfin_mesh_unitmeshes.la

libdolfin_mesh_la_SOURCES = BoundaryMesh.cpp\
                            Box.cpp\
                            Connectivity.cpp\
                            DMesh.cpp\
                            EuclideanBasis.cpp\
                            GeometricPredicates.cpp\
                            GTSInterface.cpp\
                            IntersectionDetector.cpp\
                            LoadBalancer.cpp\
                            MappedManifold.cpp\
                            Mesh.cpp\
                            MeshDependent.cpp\
                            MeshDistributedData.cpp\
                            MeshEditor.cpp\
                            MeshGeometry.cpp\
                            MeshOrdering.cpp\
                            MeshPartition.cpp\
                            MeshQuality.cpp\
                            MeshRenumber.cpp\
                            MeshSmoother.cpp\
                            MeshSmoothing.cpp\
                            MeshTopology.cpp\
                            MetisInterface.cpp\
                            MPIMeshCommunicator.cpp\
                            Rectangle.cpp\
                            RefinementPattern.cpp\
                            RivaraRefinement.cpp\
                            StructuredGrid.cpp\
                            SubDomain.cpp\
                            UniformRefinement.cpp\
                            VertexNormal.cpp\
                            ZoltanInterface.cpp
