load(sofa/pre)
defineAsPlugin(ldidetection)

DEFINES -= SOFA_GPU_CUDA

TEMPLATE=lib
TARGET = ldidetection

DEFINES += SOFA_BUILD_LDIDETECTION

SOURCES = \
         lib/DepthPeeling.cpp \
         lib/DepthPeelingUtility.cpp \
         lib/LDIDetection.cpp \
         lib/LDIPenalityContact.cpp \
         lib/LDIPenalityContactForceField.cpp \
         lib/LayeredDepthImagesPipeline.cpp \
         lib/initldidetection.cpp \
         lib/ContactConstraint.cpp \
         lib/LDIConstraintContact.cpp 

HEADERS = \
         lib/DepthPeeling.h \
         lib/DepthPeelingUtility.h \
         lib/LDIDetection.h \
         lib/LDIPenalityContact.h \
         lib/LDIPenalityContactForceField.h \
         lib/LayeredDepthImagesPipeline.h \
         lib/initldidetection.h \
         lib/ContactConstraint.h \
         lib/ContactConstraint.inl \
         lib/LDIConstraintContact.h \
         lib/LDIConstraintContact.inl \
         lib/LDIPenalityContact.inl \
         lib/LDIPenalityContactForceField.inl 

README_FILE = LICENCE.txt

#TODO: add an install target for README files

unix : QMAKE_POST_LINK = cp $$SRC_DIR/$$README_FILE $$LIB_DESTDIR 
win32 : QMAKE_POST_LINK = copy \"$$toWindowsPath($$SRC_DIR/$$README_FILE)\" \"$$LIB_DESTDIR"

load(sofa/post)
