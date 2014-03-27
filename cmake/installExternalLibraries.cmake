# -------------------------------------------------------------------------
# Find and install SOFA external libs
# -------------------------------------------------------------------------
if (WIN32)

  # For windows, just grab everything in the dependency dir and install
  file(GLOB files ${SOFA_LIB_OS_DIR}/*.dll)
  foreach(f ${files})
    install(FILES ${f} DESTINATION bin COMPONENT Runtime)
  endforeach()

else()

  # Assume that all the external libs are in the form lib_LIBRARIES and lib_INCLUDE_DIR
  set(libraries "GLUT;GLEW;PNG;ZLIB;")
  foreach(lib ${libraries})
    foreach(target ${${lib}_LIBRARIES})
      if (EXISTS ${target} AND NOT IS_DIRECTORY ${target})

        # Note: Library symlinks are always named "lib.*.dylib" on mac or "lib.so.*" on linux
        get_filename_component(lib_dir ${target} PATH)
        get_filename_component(lib_real_path ${target} REALPATH)
        get_filename_component(lib_real_dir ${lib_real_path} PATH)
        get_filename_component(lib_name ${target} NAME_WE)
        install(DIRECTORY ${lib_real_dir}/
          DESTINATION ${SOFA-INSTALL_BIN_DIR} COMPONENT Runtime
          FILES_MATCHING PATTERN "${lib_name}*"
          PATTERN "${lib_name}*.a" EXCLUDE)
        if (NOT "${lib_dir}" STREQUAL "${lib_real_dir}")
          install(CODE "message(\"External library ${lib} at ${target} sym-links to a different directory: ${lib_real_dir} .\")"
            CONFIGURATIONS Release)
          message(WARNING "External library ${lib} at ${target} sym-links to a different directory: ${lib_real_dir} .")
        endif()
      endif()
    endforeach()
  endforeach()

endif()
