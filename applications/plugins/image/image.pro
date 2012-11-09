load(sofa/pre)
defineAsPlugin(image)

TEMPLATE = lib
TARGET = sofaimage

DEFINES += SOFA_BUILD_IMAGE

# Make sure there are no cross-dependencies
INCLUDEPATH -= $$SOFA_INSTALL_INC_DIR/applications

INCLUDEPATH += $$SOFA_INSTALL_INC_DIR/extlibs

contains(DEFINES, SOFA_IMAGE_HAVE_OPENCV) { # should be "SOFA_HAVE_OPENCV" -> use "SOFA_IMAGE_HAVE_OPENCV" until the opencv plugin is fixed..
	INCLUDEPATH += $$SOFA_OPENCV_PATH
        LIBS += -lml  -lcvaux -lhighgui -lcv -lcxcore
        }

contains(DEFINES, SOFA_HAVE_LIBFREENECT) {
        INCLUDEPATH += $$SOFA_LIBFREENECT_PATH
        LIBS += -lfreenect -lfreenect_sync
        HEADERS += Kinect.h
        SOURCES += Kinect.cpp
        }

HEADERS += \
	initImage.h \
	ImageTypes.h \
	ImageContainer.h \
   	ImageViewer.h \
	ImageFilter.h \
        TransferFunction.h \
        ImageValuesFromPositions.h \
        MergeImages.h \
        ImageAccumulator.h \
        DepthMapToMeshEngine.h \
        MeshToImageEngine.h \
        MarchingCubesEngine.h \
        ImageSampler.h \
        ImageExporter.h \
	ImagePlaneWidget.h \
	ImageTransformWidget.h \
	HistogramWidget.h \
#	QImageMouseButtonsWidget.h \
	VectorVisualizationWidget.h \
        VectorVis.h \
        ImageAlgorithms.h

SOURCES += \
	initImage.cpp \
	ImageContainer.cpp \
	ImageViewer.cpp \
	ImageFilter.cpp \
        TransferFunction.cpp \
        ImageValuesFromPositions.cpp \
        MergeImages.cpp \
        ImageAccumulator.cpp \
        DepthMapToMeshEngine.cpp \
        MeshToImageEngine.cpp \
        MarchingCubesEngine.cpp \
        ImageSampler.cpp \
        ImageExporter.cpp \
	ImagePlaneWidget.cpp \
	ImageTransformWidget.cpp \
	HistogramWidget.cpp \
#	QImageMouseButtonsWidget.cpp \
	VectorVisualizationWidget.cpp \

README_FILE = image.txt
	
unix : QMAKE_POST_LINK = cp $$SRC_DIR/$$README_FILE $$LIB_DESTDIR
win32 : QMAKE_POST_LINK = copy \"$$toWindowsPath($$SRC_DIR/$$README_FILE)\" \"$$LIB_DESTDIR\"

load(sofa/post)
	
