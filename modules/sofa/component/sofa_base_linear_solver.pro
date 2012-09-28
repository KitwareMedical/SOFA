load(sofa/pre)

TEMPLATE = lib
TARGET = sofa_base_linear_solver

DEFINES += SOFA_BUILD_BASE_LINEAR_SOLVER

win32 {
	QMAKE_CXXFLAGS += /bigobj
}

HEADERS += initBaseLinearSolver.h \
           linearsolver/CGLinearSolver.h \
           linearsolver/CGLinearSolver.inl \
           linearsolver/CholeskySolver.h \
           linearsolver/CholeskySolver.inl \
           linearsolver/BTDLinearSolver.h \
           linearsolver/BTDLinearSolver.inl \
           linearsolver/FullVector.h \
           linearsolver/FullMatrix.h \
           linearsolver/DiagonalMatrix.h \
           linearsolver/SingleMatrixAccessor.h \
           linearsolver/SparseMatrix.h \
           linearsolver/CompressedRowSparseMatrix.h \
           linearsolver/CompressedRowSparseMatrix.inl \
           linearsolver/GraphScatteredTypes.h \
           linearsolver/DefaultMultiMatrixAccessor.h \
           linearsolver/MatrixLinearSolver.h \
           linearsolver/MatrixExpr.h \
           linearsolver/matrix_bloc_traits.h \
           linearsolver/BlocMatrixWriter.h

SOURCES += initBaseLinearSolver.cpp \
           linearsolver/CGLinearSolver.cpp \
           linearsolver/CholeskySolver.cpp \
           linearsolver/BTDLinearSolver.cpp \
           linearsolver/FullVector.cpp \
           linearsolver/GraphScatteredTypes.cpp \
           linearsolver/DefaultMultiMatrixAccessor.cpp \
           linearsolver/MatrixLinearSolver.cpp \
           linearsolver/SingleMatrixAccessor.cpp


# Make sure there are no cross-dependencies
INCLUDEPATH -= $$SOFA_INSTALL_INC_DIR/applications
DEPENDPATH -= $$SOFA_INSTALL_INC_DIR/applications

#exists(component-local.cfg): include(component-local.cfg)

load(sofa/post)
 
