
FFC = ffc -l dolfin

TESTS = main/MPI \
				mesh/MeshEditor \
				mesh/Mesh \
				mesh/Edge \
				mesh/GlobalFacetMap \
				mesh/BoundaryMesh \
				mesh/VertexNormal \
				fem/Assembler \
				fem/NodeNormal \
				fem/PeriodicBC \
				io/XMLFile \
				function/FunctionInterpolation \
				function/FunctionDecomposition \
				adaptivity/AdaptiveRefinement 

all: regen check

regen:
	@for d in $(TESTS); \
		do ( cd $$d && \
					( for i in *.ufl; \
							do if test -f $$i; then $(FFC) $$i; fi; done; ) && \
					cd - ); \
		done;

check:
	@set -e; \
	 for d in $(TESTS); \
		do ( cd $$d && \
					( if ! ($(MAKE) clean && $(MAKE)); \
							then echo "Compilation failed: $$d"; exit 1; fi; ) && \
					cd - ); \
		done;

