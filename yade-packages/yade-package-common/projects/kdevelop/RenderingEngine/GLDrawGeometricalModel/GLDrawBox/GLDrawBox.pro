# File generated by kdevelop's qmake manager. 
# ------------------------------------------- 
# Subdir relative project main directory: ./RenderingEngine/GLDrawGeometricalModel/GLDrawBox
# Target is a library:  

HEADERS += GLDrawBox.hpp 
SOURCES += GLDrawBox.cpp 
LIBS += -lBox \
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
