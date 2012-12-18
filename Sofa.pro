message( "PRE-CONFIG: " $${CONFIG})
load(sofa/pre)

TEMPLATE = subdirs

!contains(DEFINES, SOFA_DEV): message("WARNING: SOFA_DEV not defined, in-development code will be disabled!")

contains(DEFINES, SOFA_RELEASE): message("WARNING: SOFA_RELEASE defined, in-development code will be disabled!")

#message( "artifacts registry : $$artifacts_registry" )


########################################################################
# Generate SUBDIRS specifications to build everything
########################################################################

buildEnabledArtifacts()

########################################################################
# Print current config
########################################################################

message( "====== SOFA Build Configuration ======")

contains(DEFINES,SOFA_DEV){ # BEGIN SOFA_DEV
message( "==== UNSTABLE DEVELOPMENT VERSION ====")
} # END SOFA_DEV

win32 {
  message( "|  Platform: Windows")
}
else:macx {
  message( "|  Platform: MacOS")
}
else:unix {
  message( "|  Platform: Linux/Unix")
}

win32:!build_pass {
	message( "|  Mode : PROJECT")
}
else {
contains (CONFIGDEBUG, debug) {
	contains( CONFIGSTATIC, static) {
		message( "|  Mode: DEBUG with static compilation")
	}
	else {
	  message( "|  Mode: DEBUG")
	}
}
contains (CONFIGDEBUG, release) {
  contains (QMAKE_CXXFLAGS,-g) {
    message( "|  Mode: RELEASE with debug symbols")
  }
  else {
    contains (CONFIGDEBUG, profile) {
      message( "|  Mode: RELEASE with profiling")
    }
    else {
			contains (CONFIGSTATIC, static) {
	      message( "|  Mode: RELEASE with static compilation")
			}
			else {
				message( "|  Mode : RELEASE")
			}
    }
  }
}
}



contains(DEFINES,SOFA_QT4) {
  message( "|  Qt version: 4.x")
}
else {
  message( "|  Qt version: 3.x")
}

contains(DEFINES,SOFA_RDTSC) {
  message( "|  RDTSC timer: ENABLED")
}
else {
  message( "|  RDTSC timer: DISABLED")
}

contains(DEFINES,SOFA_HAVE_BOOST) {
  message( "|  BOOST libraries: ENABLED")
}
else {
  message( "|  BOOST libraries: DISABLED")
}

contains(DEFINES,SOFA_HAVE_DAG) {
  message( "|  Directed Acyclic Graph (DAG): ENABLED")
}
else {
  message( "|  Directed Acyclic Graph (DAG): DISABLED")
}

contains(DEFINES,SOFA_HAVE_BGL) {
  message( "|  Boost Graph Library (BGL): ENABLED")
}
else {
  message( "|  Boost Graph Library (BGL): DISABLED")
}

contains(DEFINES,SOFA_HAVE_PYTHON) {
  message( "|  PYTHON script support: ENABLED")
}
else {
  message( "|  PYTHON script support: DISABLED")
}

contains(DEFINES,SOFA_XML_PARSER_TINYXML) {
  message( "|  TinyXML parser: ENABLED")
}
else {
  message( "|  TinyXML parser: DISABLED")
}

contains(DEFINES,SOFA_XML_PARSER_LIBXML) {
  message( "|  LibXML parser: ENABLED")
}
else {
  message( "|  LibXML parser: DISABLED")
}

contains(DEFINES,SOFA_HAVE_PNG) {
  message( "|  PNG support: ENABLED")
}
else {
  message( "|  PNG support: DISABLED")
}

contains(DEFINES,SOFA_HAVE_GLEW) {
  message( "|  OpenGL Extensions support using GLEW: ENABLED")
}
else {
  message( "|  OpenGL Extensions support using GLEW: DISABLED")
}

contains(DEFINES,SOFA_GPU_CUDA) {
  message( "|  GPU support using CUDA: ENABLED")
}
else {
  message( "|  GPU support using CUDA: DISABLED")
}
contains(DEFINES,SOFA_SMP) {
  message( "|   Sofa-Parallel: ENABLED ")
  message( "| KAAPI_DIR=$${KAAPI_DIR}")
}
else {
  message( "|  Sofa-Parallel: DISABLED")
}

contains(DEFINES,SOFA_GPU_OPENCL) {
  message( "|  GPU support using OPENCL: ENABLED")
}
else {
  message( "|  GPU support using OPENCL: DISABLED")
}

contains(DEFINES,SOFA_PML) {
  message( "|  PML/LML support: ENABLED")
}
else {
  message( "|  PML/LML support: DISABLED")
}


contains(DEFINES,SOFA_HAVE_CSPARSE) {
  message( "|  CSPARSE library : ENABLED")
}
else {
  message( "|  CSPARSE library : DISABLED")
}

contains(DEFINES,SOFA_HAVE_METIS) {
  message( "|  METIS library : ENABLED")
}
else {
  message( "|  METIS library : DISABLED")
}

contains(DEFINES,SOFA_HAVE_TAUCS) {
  message( "|  TAUCS library : ENABLED")
contains(DEFINES,SOFA_HAVE_CILK) {
  message( "|  CILK library : ENABLED")
} else {
  message( "|  CILK library : DISABLE")
}
}
else {
  message( "|  TAUCS library : DISABLED")
}


contains(DEFINES,SOFA_GUI_GLUT) {
  message( "|  GLUT GUI: ENABLED")
}
else {
  message( "|  GLUT GUI: DISABLED")
}

!contains(DEFINES,SOFA_GUI_QTVIEWER) {
!contains(DEFINES,SOFA_GUI_QGLVIEWER) {
{
  message( "|  Qt GUI: DISABLED")
}
#else {
 # message( "|  Qt GUI: ENABLED")
#}
}
else {
  message( "|  Qt GUI: ENABLED")
}
}
else {
  message( "|  Qt GUI: ENABLED")
}

contains(DEFINES,SOFA_GUI_QTVIEWER) {
  message( "|  -  Qt OpenGL viewer: ENABLED")
}
else {
  message( "|  -  Qt OpenGL viewer: DISABLED")
}

contains(DEFINES,SOFA_GUI_QGLVIEWER) {
  message( "|  -  Qt QGLViewer viewer: ENABLED")
}
else {
  message( "|  -  Qt QGLViewer viewer: DISABLED")
}

message( "======================================")
message( "|  CONFIG: " $${CONFIG})
message( "|  DEFINES: " $${DEFINES})
message( "======================================")



unix {
  contains(DEFINES, SOFA_QT4):DOLLAR="\\$"
  !contains(DEFINES, SOFA_QT4):DOLLAR="\$"
  contains(DEFINES, SOFA_SMP) {
    system(echo "export SOFA_DIR=$${PWD}" >config-Sofa-parallel.sh)
    system(echo "export KAAPI_DIR=$${KAAPI_DIR}" >>config-Sofa-parallel.sh)
    system(echo "export LD_LIBRARY_PATH=$${DOLLAR}SOFA_DIR/lib/linux:$${DOLLAR}KAAPI_DIR/lib:$${DOLLAR}LD_LIBRARY_PATH" >>config-Sofa-parallel.sh)
    system(echo "export PATH=$${DOLLAR}SOFA_DIR/bin:$${DOLLAR}KAAPI_DIR/bin:$${DOLLAR}PATH" >>config-Sofa-parallel.sh)
    contains(DEFINES, SOFA_GPU_CUDA) {
      system(echo "export CUDA_DIR=$${CUDA_DIR}" >>config-Sofa-parallel.sh)
      system(echo "export LD_LIBRARY_PATH=$${DOLLAR}CUDA_DIR/lib:$${DOLLAR}CUDA_DIR/lib64:$${DOLLAR}LD_LIBRARY_PATH" >>config-Sofa-parallel.sh)
      system(echo "export PATH=$${DOLLAR}CUDA_DIR/bin:$${DOLLAR}PATH" >>config-Sofa-parallel.sh)
    }
  }
}


########################################################################
# Export activated libraries and their dependencies as a graphviz file.
########################################################################
unix {
#  outputBuildGraph(Sofa-build.dot)
#  message(Generating Sofa-build.pdf)
#  system(dot -Tpdf -oSofa-build.pdf Sofa-build.dot)
}

load(sofa/post)

# print all SOFA DEFINES into a standard file format (mainly used for build with CMake)
message("Write temporarily the sofa DEFINES in sofaDefines$${LIBSUFFIX}.cfg to let CMake get them")
# under windows, lets try also: copy /y NUL EmptyFile.txt >NUL OR echo. 2>EmptyFile.txt to clear the textFile
win32 { system( type NUL > sofaDefines$${LIBSUFFIX}.cfg ) }
win32 { system( for %G in ($${DEFINES}) do echo %G>>sofaDefines$${LIBSUFFIX}.cfg ) }
unix  { system( echo >sofaDefines$${LIBSUFFIX}.cfg ) }
unix  { system( for define in $${DEFINES}; do echo $define>>sofaDefines$${LIBSUFFIX}.cfg; done ) }
