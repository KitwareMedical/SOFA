include(${SOFA_CMAKE_DIR}/preProject.cmake)

set(HEADER_FILES

	lib/DepthPeeling.h
	lib/DepthPeelingUtility.h
	lib/LDIDetection.h
	lib/LDIPenalityContact.h
	lib/LDIPenalityContactForceField.h
	lib/LayeredDepthImagesPipeline.h
	lib/initldidetection.h
	lib/ContactConstraint.h
	lib/ContactConstraint.inl
	lib/LDIConstraintContact.h
	lib/LDIConstraintContact.inl
	lib/LDIPenalityContact.inl
	lib/LDIPenalityContactForceField.inl 
	)

set(SOURCE_FILES

    lib/DepthPeeling.cpp
	lib/DepthPeelingUtility.cpp
	lib/LDIDetection.cpp
	lib/LDIPenalityContact.cpp
	lib/LDIPenalityContactForceField.cpp
	lib/LayeredDepthImagesPipeline.cpp
	lib/initldidetection.cpp
	lib/ContactConstraint.cpp
	lib/LDIConstraintContact.cpp
	)
	
set(README_FILES

    LICENCE.txt
	)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

AddCompilerDefinitions("SOFA_BUILD_LDIDETECTION" )
AddLinkerDependencies(SofaGuiQt)


include(${SOFA_CMAKE_DIR}/postProject.cmake)
