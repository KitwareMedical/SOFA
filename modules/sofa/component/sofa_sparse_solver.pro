load(sofa/pre)

TEMPLATE = lib
TARGET = sofa_sparse_solver

DEFINES += SOFA_BUILD_SPARSE_SOLVER

HEADERS += initSparseSolver.h \
           linearsolver/PrecomputedLinearSolver.h \
           linearsolver/PrecomputedLinearSolver.inl \
           linearsolver/SparseCholeskySolver.h \
           linearsolver/SparseLUSolver.h \
           linearsolver/SparseLDLSolverImpl.h \
           linearsolver/SparseLDLSolver.h \
           linearsolver/SparseLDLSolver.inl

SOURCES += initSparseSolver.cpp \
           linearsolver/PrecomputedLinearSolver.cpp \
           linearsolver/SparseCholeskySolver.cpp \
           linearsolver/SparseLUSolver.cpp \
           linearsolver/SparseLDLSolver.cpp


# Make sure there are no cross-dependencies
INCLUDEPATH -= $$SOFA_INSTALL_INC_DIR/applications
DEPENDPATH -= $$SOFA_INSTALL_INC_DIR/applications

contains(DEFINES, SOFA_HAVE_METIS) : !contains(DEFINES, SOFA_EXTLIBS_METIS) {
LIBS += -lmetis
}

#exists(component-local.cfg): include(component-local.cfg)

load(sofa/post)
 
