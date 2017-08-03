# Written by Emanuele Rodola
# USI Lugano, Apr 2016
# -Wl,--out-implib,./$(LIBNAME).a

LIBNAME=geodesic_matlab_api.so
BOOST_DIR=/usr/include/boost/
CXXFLAGS=-Wall -O2 -fPIC

all:
	g++ -c ./geodesic_matlab_api.cpp -o ./geodesic_matlab_api.o $(CXXFLAGS) -I$(BOOST_DIR)
	g++ -shared -o ./$(LIBNAME) ./geodesic_matlab_api.o

clean:
	rm -f ./geodesic_matlab_api.o ./$(LIBNAME) ./$(LIBNAME).a
