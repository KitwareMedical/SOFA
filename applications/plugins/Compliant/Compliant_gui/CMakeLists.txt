include(${SOFA_CMAKE_DIR}/preProject.cmake)



set(HEADER_FILES

	initCompliant_gui.h
	CompliantAttachPerformer.h
	CompliantAttachPerformer.inl
	
	)
	
	
set(SOURCE_FILES

	initCompliant_gui.cpp
	CompliantAttachPerformer.cpp

	)


if(WIN32)
	include_directories("${SOFA_EXTLIBS_DIR}/SuiteSparse/cholmod/Include")
	link_directories("${SOFA_EXTLIBS_DIR}/SuiteSparse/cholmod/Lib")
else()
	include_directories("/usr/include/suitesparse/")
endif()

	
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})

AddCompilerDefinitions("SOFA_BUILD_COMPLIANT_GUI")
AddLinkerDependencies(SofaGuiMain Compliant)


include(${SOFA_CMAKE_DIR}/postProject.cmake)
