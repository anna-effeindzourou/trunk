# File generated by kdevelop's qmake manager. 
# ------------------------------------------- 
# Subdir relative project main directory: ./RenderingEngine/GLDrawGeometricalModel/GLDrawMesh2D
# Target is a library:  

HEADERS += GLDrawMesh2D.hpp 
SOURCES += GLDrawMesh2D.cpp 
LIBS += -lMesh2D \
        -lyade-lib-opengl \
        -rdynamic 
QMAKE_LIBDIR = ../../../../../bin \
               /usr/local/lib/yade/yade-libs/ 
QMAKE_CXXFLAGS_RELEASE += -lpthread \
                          -pthread 
QMAKE_CXXFLAGS_DEBUG += -lpthread \
                        -pthread 
DESTDIR = ../../../../../bin 
CONFIG += debug \
          warn_on \
          dll 
TEMPLATE = lib 
