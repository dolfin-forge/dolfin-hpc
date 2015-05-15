CFLAGS  = `pkg-config --cflags dolfin`
CFLAGS += -I../../include
CFLAGS += -g -O3
LIBS    = `pkg-config --libs dolfin`
CXX     = `pkg-config --variable=compiler dolfin`

DEST    = demo
OBJECTS = main.o
VTKFILES= *.pvd *.vtu *.pvtu

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(VTKFILES) $(OBJECTS) $(DEST)

$(DEST): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -c $<
