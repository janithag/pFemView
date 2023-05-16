# Makefile - author Janitha Gunatilake
#************************************************************
#
#
#include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
#
#
OBJS =  pfemvtk.o element.o hex.o hexh.o wedge.o wedgeh.o tet.o teth.o quad.o quadh.o tri.o trih.o line.o lineh.o hFun.o trihFun.o triFun.o tetraFun.o linFun.o linehFun.o hVisualize.o
CC = g++ 
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

pfemvtk.exe: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o pfemvtk.exe

pfemvtk.o: hVisualize.hpp pfemvtk.cpp
	$(CC) $(CFLAGS) pfemvtk.cpp

hVisualize.o: hVisualize.hpp hVisualize.cpp
	$(CC) $(CFLAGS) hVisualize.cpp

element.o: element.hpp element.cpp basis.hpp basish.hpp
	$(CC) $(CFLAGS) element.cpp

hex.o: basis.hpp hex.cpp friendFuns.hpp
	$(CC) $(CFLAGS) hex.cpp

hexh.o: basish.hpp hexh.cpp main.hpp
	$(CC) $(CFLAGS) hexh.cpp

wedge.o: basis.hpp wedge.cpp friendFuns.hpp
	$(CC) $(CFLAGS) wedge.cpp

wedgeh.o: basish.hpp wedgeh.cpp main.hpp
	$(CC) $(CFLAGS) wedgeh.cpp

tet.o: basis.hpp tet.cpp friendFuns.hpp
	$(CC) $(CFLAGS) tet.cpp

teth.o: basish.hpp teth.cpp
	$(CC) $(CFLAGS) teth.cpp      

quad.o: basis.hpp quad.cpp friendFuns.hpp
	$(CC) $(CFLAGS) quad.cpp

quadh.o: basish.hpp quadh.cpp
	$(CC) $(CFLAGS) quadh.cpp

tri.o: basis.hpp tri.cpp friendFuns.hpp
	$(CC) $(CFLAGS) tri.cpp

trih.o: basish.hpp trih.cpp
	$(CC) $(CFLAGS) trih.cpp

line.o: basis.hpp line.cpp friendFuns.hpp 
	$(CC) $(CFLAGS) line.cpp

lineh.o: basish.hpp lineh.cpp
	$(CC) $(CFLAGS) lineh.cpp

hFun.o: hFun.cpp
	$(CC) $(CFLAGS) hFun.cpp

trihFun.o: hFun.hpp trihFun.cpp
	$(CC) $(CFLAGS) trihFun.cpp

triFun.o: triFun.cpp
	$(CC) $(CFLAGS) triFun.cpp

tetraFun.o: tetraFun.cpp friendFuns.hpp
	$(CC) $(CFLAGS) tetraFun.cpp

linFun.o: linFun.cpp
	$(CC) $(CFLAGS) linFun.cpp

linehFun.o: hFun.hpp linehFun.cpp
	$(CC) $(CFLAGS) linehFun.cpp

clean:
	rm -f $(OBJS) pfemvtk.exe

