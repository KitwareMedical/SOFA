load(sofa/pre)

TEMPLATE = lib
TARGET = sofa_opengl_visual

DEFINES += SOFA_BUILD_OPENGL_VISUAL

HEADERS += initOpenGLVisual.h \
           visualmodel/OglModel.h \
           visualmodel/OglViewport.h \
           visualmodel/Light.h \
           visualmodel/LightManager.h \
           visualmodel/PointSplatModel.h \
           visualmodel/OglRenderingSRGB.h \
           visualmodel/ClipPlane.h \
		   visualmodel/CompositingVisualLoop.h

SOURCES += initOpenGLVisual.cpp \
           visualmodel/OglModel.cpp \
           visualmodel/OglViewport.cpp \
           visualmodel/Light.cpp \
           visualmodel/LightManager.cpp \
           visualmodel/PointSplatModel.cpp \
           visualmodel/OglRenderingSRGB.cpp \
           visualmodel/ClipPlane.cpp \
           visualmodel/CompositingVisualLoop.cpp

contains(DEFINES,SOFA_HAVE_GLEW) {    # seems to be built without GLEW define since we include glext.h from helper/system/gl.h
HEADERS += visualmodel/OglAttribute.h \
           visualmodel/OglAttribute.inl \
           visualmodel/OglShader.h \
           visualmodel/OglShaderMacro.h \
           visualmodel/OglShaderVisualModel.h \
           visualmodel/OglShadowShader.h \
           visualmodel/OglTetrahedralModel.h \
           visualmodel/OglTetrahedralModel.inl \
           visualmodel/OglTexture.h \
           visualmodel/OglVariable.h \
           visualmodel/OglVariable.inl \
           visualmodel/PostProcessManager.h \
           visualmodel/SlicedVolumetricModel.h \
   		   visualmodel/VisualManagerPass.h \
		   visualmodel/VisualManagerSecondaryPass.h

SOURCES += visualmodel/OglShader.cpp \
           visualmodel/OglShaderMacro.cpp \
           visualmodel/OglShaderVisualModel.cpp \
           visualmodel/OglShadowShader.cpp \
           visualmodel/OglTetrahedralModel.cpp \
           visualmodel/OglTexture.cpp \
           visualmodel/OglVariable.cpp \
           visualmodel/PostProcessManager.cpp \
           visualmodel/SlicedVolumetricModel.cpp \
		   visualmodel/VisualManagerPass.cpp  \
		   visualmodel/VisualManagerSecondaryPass.cpp
}

# Make sure there are no cross-dependencies
INCLUDEPATH -= $$SOFA_INSTALL_INC_DIR/applications
DEPENDPATH -= $$SOFA_INSTALL_INC_DIR/applications

#exists(component-local.cfg): include(component-local.cfg)

load(sofa/post)
 
