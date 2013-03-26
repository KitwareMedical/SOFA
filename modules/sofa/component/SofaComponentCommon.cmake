cmake_minimum_required(VERSION 2.8)

project("SofaComponentCommon")

include(${SOFA_CMAKE_DIR}/pre.cmake)

set(HEADER_FILES

    initComponentCommon.h

    )
    
set(SOURCE_FILES

    initComponentCommon.cpp

    )
    
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})

set(COMPILER_DEFINES "SOFA_BUILD_COMPONENT_COMMON")
set(LINKER_DEPENDENCIES SofaLoader SofaDeformable SofaSimpleFem SofaObjectInteraction SofaExplicitOdeSolver SofaImplicitOdeSolver SofaEigen2Solver SofaMeshCollision)

include(${SOFA_CMAKE_DIR}/post.cmake)