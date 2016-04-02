CFLAGS  = `pkg-config --cflags dolfin`
LIBS    = `pkg-config --libs dolfin`
CXX     = `pkg-config --variable=compiler dolfin`
NPROCS ?= 1
TARGET  = exe
OBJECTS = main.o

all: $(TARGET)

install:

clean:
	-rm -f *.o core *.core $(OBJECTS) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -c $<
	
test:	clean $(TARGET)
	mpirun -n $(NPROCS) ./$(TARGET)
