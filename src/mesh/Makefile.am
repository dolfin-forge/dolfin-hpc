SUBDIRS=celltypes entities mappings unitmeshes utilities

AM_CPPFLAGS = -I@top_srcdir@/include
noinst_LTLIBRARIES = libdolfin_mesh.la

libdolfin_mesh_la_LIBADD = celltypes/libdolfin_mesh_celltypes.la\
                           entities/libdolfin_mesh_entities.la\
                           mappings/libdolfin_mesh_mappings.la\
                           unitmeshes/libdolfin_mesh_unitmeshes.la\
                           utilities/libdolfin_mesh_utilities.la

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
                            MeshTopology.cpp\
                            MetisInterface.cpp\
                            Rectangle.cpp\
                            RefinementPattern.cpp\
                            StructuredGrid.cpp\
                            SubDomain.cpp\
                            VertexNormal.cpp\
                            ZoltanInterface.cpp
