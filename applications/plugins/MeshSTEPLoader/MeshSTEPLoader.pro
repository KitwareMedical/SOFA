######  PLUGIN TARGET
load(sofa/pre)
defineAsPlugin(MeshSTEPLoader)

TARGET = MeshSTEPLoader

###### SPECIFIC PLUGIN CONFIGURATION, you should modify it to configure your plugin

DEFINES += SOFA_BUILD_MESHSTEPLOADERPLUGIN

SOURCES = \
MeshSTEPLoader.cpp \
ParametricTriangleTopologyContainer.cpp \
SingleComponent.cpp \
STEPShapeMapping.cpp \
initMeshSTEPLoader.cpp

HEADERS = \
MeshSTEPLoader.h\
ParametricTriangleTopologyContainer.h \
SingleComponent.inl\
SingleComponent.h \
STEPShapeMapping.h \
initMeshSTEPLoader.h



README_FILE = PluginMeshSTEPLoader.txt

unix {
INCLUDEPATH += /usr/include/opencascade  # FF: is this really useful ?
DEPENDPATH += /usr/include/opencascade
LIBS += -lTKernel -lTKMath -lTKAdvTools -lGL -lTKG2d -lTKG3d -lTKGeomBase -lTKBRep -lTKGeomAlgo -lTKTopAlgo -lTKPrim -lTKBO -lTKHLR -lTKMesh -lTKShHealing -lTKBool -lTKXMesh -lTKFillet -lTKFeat -lTKOffset -lTKSTL -lTKXSBase -lTKSTEPBase -lTKIGES -lTKSTEPAttr -lTKSTEP209 -lTKSTEP    -lTKService -lTKV2d -lTKV3d -lTKOpenGl -lTKMeshVS -lTKNIS -lTKVRML
}

win32 {
INCLUDEPATH += $$OPEN_CASCADE_DIR/inc
DEPENDPATH += $$OPEN_CASCADE_DIR/inc
LIBS += -l$$OPEN_CASCADE_DIR/win32/lib/*
QMAKE_CXXFLAGS += /DWNT
}
unix : QMAKE_POST_LINK = cp $$SRC_DIR/$$README_FILE $$LIB_DESTDIR 
win32 : QMAKE_POST_LINK = copy \"$$toWindowsPath($$SRC_DIR/$$README_FILE)\" \"$$LIB_DESTDIR\"

load(sofa/post)
