include_directories(${PROJECT_BINARY_DIR}/include)

if (ECOSIM_NEED_PFUNIT EQUAL 1)
# Build pfunit, a fortran unit testing library.
  if (NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libpfunit.a)
    message("Configuring pfunit...")
    execute_process(COMMAND cmake -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/pfunit
      RESULT_VARIABLE stat)

    if (HAVE_MPI EQUAL 1)
      message(FATAL_ERROR "BeTR pUFnit MPI build is not yet supported.")

      set(PFUNIT_CMAKE_OPTS ${PFUNIT_CMAKE_OPTS} -DMPI=YES -DMPI_MPIFC=${CMAKE_Fortran_COMPILER})
    else()
      set(PFUNIT_CMAKE_OPTS ${PFUNIT_CMAKE_OPTS} -DMPI=NO)
    endif()

    set(PFUNIT_CMAKE_OPTS
      -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    )
    execute_process(COMMAND cmake
      ${PFUNIT_CMAKE_OPTS}
      ${CMAKE_CURRENT_SOURCE_DIR}/pfunit
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pfunit
      OUTPUT_VARIABLE pfunit_config.log ERROR_VARIABLE pfunit_config_err
      RESULT_VARIABLE stat)

    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/pfunit_config.log ${pfunit_config_log})

    if (NOT stat EQUAL 0)
      if (TRAVIS_CI)
        message("Error configuring pfunit:")
        message(${pfunit_config_log})
        message(FATAL_ERROR ${pfunit_config_err})
      else()
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/pfunit_config_errors.log ${pfunit_config_err})
        message(FATAL_ERROR "Configuration of pfunit library failed.")
      endif()
    endif()

    message("Building and installing pfunit...")
    execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pfunit
                  OUTPUT_VARIABLE pfunit_build_log ERROR_VARIABLE pfunit_build_err
                  RESULT_VARIABLE stat)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/pfunit_build.log ${pfunit_build_log})
    if (NOT stat EQUAL 0)
      if (TRAVIS_CI)
        message("Error building pfunit:")
        message(${pfunit_build_log})
        message(FATAL_ERROR ${pfunit_build_err})
      else()
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/pfunit_build_errors.log ${pfunit_build_err})
        message(FATAL_ERROR "Build of pfunit library failed.")
      endif()
    endif()
  endif()
  add_library(pfunit STATIC IMPORTED GLOBAL)
  set_target_properties(pfunit PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libpfunit.a)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/3rd-party/pfunit/include/" PARENT_SCOPE)
  set(ENV{PFUNIT} ${PROJECT_BINARY_DIR} CACHE)
endif()


#remove netcdf dependence
if (ECOSIM_HAVE_NETCDF EQUAL 1)
return ()
endif()

# Build zlib for compression.
if (NOT EXISTS ${Z_LIBRARY})
  set(ZLIB_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --static)

  message("Preparing zlib...")
  message("PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")
  message("CMAKE_COMMAND: ${CMAKE_COMMAND}")
  message("CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
  message("CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")

  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/zlib ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Configuring zlib...")
  message("CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
  message("CMAKE_C_FLAGS: ${CMAKE_C_FLAGS}")
  message("ZLIB_CONFIG_OPTS: ${ZLIB_CONFIG_OPTS}") 

  execute_process(COMMAND env CC=${CMAKE_C_COMPILER} env CFLAGS=${CMAKE_C_FLAGS} ./configure ${ZLIB_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  OUTPUT_VARIABLE zlib_config_log ERROR_VARIABLE zlib_config_err
                  RESULT_VARIABLE stat)

  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_config.log ${zlib_config_log})

  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Configuration of zlib library failed:")
      message(${zlib_config_log})
      message(FATAL_ERROR ${zlib_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_config_errors.log ${zlib_config_err})
      message(FATAL_ERROR "Configuration of zlib library failed. See ${CMAKE_CURRENT_BINARY_DIR}/zlib_config_*.log.")
    endif()
  endif()

  message("Building and installing zlib (${NUM_BUILD_THREADS} threads)...")

  execute_process(COMMAND make install -j${NUM_BUILD_THREADS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/zlib
                  OUTPUT_VARIABLE zlib_build_log ERROR_VARIABLE zlib_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_build.log ${zlib_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Build of zlib library failed:")
      message(${zlib_build_log})
      message(FATAL_ERROR ${zlib_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/zlib_build_errors.log ${zlib_build_err})
      message(FATAL_ERROR "Build of zlib library failed. See ${CMAKE_CURRENT_BINARY_DIR}/zlib_build_*.log.")
    endif()
  endif()
  message(STATUS "Z_INCLUDE_DIR ${Z_INCLUDE_DIR}")
  message(STATUS "Z_LIBRARY ${Z_LIBRARY}")
endif()
set(ECOSIM_TP_LIBS z;${ECOSIM_TP_LIBS})


# Build the curl library
if (EXISTS ${CURL_LIBRARY})
  message("Preparing curl...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/curl ${CMAKE_CURRENT_BINARY_DIR}/curl
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()
  message("Configuring curl library...")
  execute_process(COMMAND autoreconf -if WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/curl)

  execute_process(COMMAND env CC=${CMAKE_C_COMPILER} env CFLAGS=${CMAKE_C_FLAGS} ./configure --prefix ${PROJECT_BINARY_DIR} --without-ssl --without-zlib
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/curl
                  OUTPUT_VARIABLE curl_config_log ERROR_VARIABLE curl_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/curl_config.log ${curl_config_log})

  message("Building and installing curl library...")
  execute_process(COMMAND make
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/curl
                  OUTPUT_VARIABLE curl_build_log ERROR_VARIABLE curl_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/curl_build.log ${curl_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Error building curl:")
      message(${curl_build_log})
      message(FATAL_ERROR ${curl_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/curl_build_errors.log ${curl_build_err})
      message(FATAL_ERROR "Error building curl library. See ${CMAKE_CURRENT_BINARY_DIR}/curl_build*.log for details.")
    endif()
  endif()

  execute_process(COMMAND make install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/curl
                  OUTPUT_VARIABLE curl_build_log ERROR_VARIABLE curl_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/curl_build.log ${curl_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Error building curl:")
      message(${curl_build_log})
      message(FATAL_ERROR ${curl_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/curl_build_errors.log ${curl_build_err})
      message(FATAL_ERROR "Error building curl library. See ${CMAKE_CURRENT_BINARY_DIR}/curl_build*.log for details.")
    endif()
  endif()
  message(STATUS "CURL_INCLUDE_DIR ${CURL_INCLUDE_DIR}")
  message(STATUS "CURL_LIBRARY ${CURL_LIBRARY}")

  set(ECOSIM_TP_LIBS curl;${ECOSIM_TP_LIBS})
endif()

if (EXISTS ${BLOSC_LIBRARY})
  message("Preparing curl...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/c-blosc ${CMAKE_CURRENT_BINARY_DIR}/c-blosc
                  RESULT_VARIABLE stat)

  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Configuring blosc library...")
  execute_process(COMMAND cmake -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/c-blosc
                  OUTPUT_VARIABLE blosc_build_log ERROR_VARIABLE blosc_build_err
                  RESULT_VARIABLE stat)

  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/blosc_config.log ${blosc_config_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Error configuring blosc:")
      message(${blosc_config_log})
      message(FATAL_ERROR ${blosc_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/blosc_config_errors.log ${blosc_config_err})
      message(FATAL_ERROR "Error configuring blosc library. See ${CMAKE_CURRENT_BINARY_DIR}/blosc_config*_log for details.")
    endif()
  endif()
  message("Building blosc library...")
  execute_process(COMMAND cmake --build . --target install
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/c-blosc
                  OUTPUT_VARIABLE blosc_build_log ERROR_VARIABLE blosc_build_err
                  RESULT_VARIABLE stat)

  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/blosc_build_log ${blosc_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Error building blosc:")
      message(${blosc_build_log})
      message(FATAL_ERROR ${blosc_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/blosc_build_errors.log ${blosc_build_err})
      message(FATAL_ERROR "Error building blosc library. See ${CMAKE_CURRENT_BINARY_DIR}/blosc_build*_log for details.")
    endif()
  endif()
  message(STATUS "BLOSC_INCLUDE_DIR ${BLOSC_INCLUDE_DIR}")
  message(STATUS "BLOSC_LIBRARY ${BLOSC_LIBRARY}")
endif()

if (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/hdf5)
  message("Preparing hdf5...")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/hdf5 ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()
  message("Configuring hdf5...")

  if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/hdf5)
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/hdf5/CMakeCache.txt)
  else()
    execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/hdf5)
  endif()
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(HDF5_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --with-zlib=${PROJECT_BINARY_DIR} --disable-shared --enable-hl --enable-build-mode=debug)
  else()
    set(HDF5_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR} --with-zlib=${PROJECT_BINARY_DIR} --disable-shared --enable-hl)
  endif()

  execute_process(COMMAND autoreconf -if WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hdf5)

  execute_process(COMMAND env CC=${CMAKE_C_COMPILER} env CFLAGS=${CMAKE_C_FLAGS} env FC=${CMAKE_Fortran_COMPILER} env CXX=${CMAKE_CXX_COMPILER}  ./configure ${HDF5_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hdf5
                  OUTPUT_VARIABLE hdf5_config_log ERROR_VARIABLE hdf5_config_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config.log ${hdf5_config_log})

  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config.log ${hdf5_config_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Configuration of hdf5 library failed:")
      message(${hdf5_config_log})
      message(FATAL_ERROR ${hdf5_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config_errors.log ${hdf5_config_err})
      message(FATAL_ERROR "Configuration of hdf5 library failed. See ${CMAKE_CURRENT_BINARY_DIR}/hdf5_config_*.log.")
    endif()
  endif()

  message("Building and installing hdf5 (${NUM_BUILD_THREADS} threads)...")
  execute_process(COMMAND make install -j${NUM_BUILD_THREADS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/hdf5
    OUTPUT_VARIABLE hdf5_build_log ERROR_VARIABLE hdf5_build_err
    RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_build.log ${hdf5_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Build of hdf5 library failed:")
      message(${hdf5_build_log})
      message(FATAL_ERROR ${hdf5_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/hdf5_build_errors.log ${hdf5_build_err})
      message(FATAL_ERROR "Build of hdf5 library failed. See ${CMAKE_CURRENT_BINARY_DIR}/hdf5_build_*.log.")
    endif()
  endif()

  # On Macs we need to change the rpath of shared libraries so we can install them elsewhere.
  message(STATUS "${APPLE} and ${BUILD_SHARED_LIBS}")
  if (APPLE AND BUILD_SHARED_LIBS)
    execute_process(COMMAND otool -DX ${PROJECT_BINARY_DIR}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}
      OUTPUT_VARIABLE libhdf5_id)
    string(STRIP ${libhdf5_id} libhdf5_id)
    execute_process(COMMAND install_name_tool -id @rpath/lib${HDF5_LIB_NAME}${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/${HDF5_LIB_NAME}${LIB_SUFFIX})
    execute_process(COMMAND install_name_tool -id @rpath/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX})
    execute_process(COMMAND install_name_tool -change ${libhdf5_id} @rpath/libhdf5${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/${HDF5_HL_LIB_NAME}${LIB_SUFFIX})
  endif()

endif()

message(STATUS "${PROJECT_BINARY_DIR}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}")
if (HDF5_LIBRARY STREQUAL "${PROJECT_BINARY_DIR}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}")
  message(STATUS "set ECOSIM_hdf5_library from 3rd-party ${CMAKE_INSTALL_PREFIX}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}")
  set(ECOSIM_HDF5_LIBRARY "${CMAKE_INSTALL_PREFIX}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}" PARENT_SCOPE)
  set(ECOSIM_HDF5_LIBRARY "${PROJECT_BINARY_DIR}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}")

  set(ECOSIM_HDF5_HL_LIBRARY "${CMAKE_INSTALL_PREFIX}/lib/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX}" PARENT_SCOPE)
  set(ECOSIM_HDF5_HL_LIBRARY "${PROJECT_BINARY_DIR}/lib/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX}")
  set(ECOSIM_HDF5_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include" PARENT_SCOPE)
  set(ECOSIM_HDF5_INCLUDE_DIR "${PROJECT_BINARY_DIR}/include")
else()
  message(STATUS "set ECOSIM_hdf5_library from system")
  set(ECOSIM_HDF5_LIBRARY ${HDF5_LIBRARY} PARENT_SCOPE)
  set(ECOSIM_HDF5_HL_LIBRARY ${HDF5_HL_LIBRARY} PARENT_SCOPE)
  set(ECOSIM_HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIR} PARENT_SCOPE)
endif()

# If our HDF5 library is named something weird, copy it to a file with the
# regular boring name. This prevents libraries like Silo and Exodus (which
# take a while to catch up to recent developments) from freaking out.
if (NOT HDF5_LIB_NAME STREQUAL "hdf5" AND HDF5_LIBRARY MATCHES "${PROJECT_BINARY_DIR}/lib" AND NOT EXISTS ${PROJECT_BINARY_DIR}/lib/libhdf5${LIB_SUFFIX})
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${HDF5_LIBRARY} ${PROJECT_BINARY_DIR}/lib/libhdf5${LIB_SUFFIX})
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${HDF5_HL_LIBRARY} ${PROJECT_BINARY_DIR}/lib/libhdf5_hl${LIB_SUFFIX})
endif()

set(ECOSIM_TP_LIBS ${HDF5_LIBRARIES};${ECOSIM_TP_LIBS})
if (LINUX EQUAL 1)
  # Needed to link hdf5 / netcdf
  set(ECOSIM_TP_LIBS "-ldl";${ECOSIM_TP_LIBS})
  #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ldl" PARENT_SCOPE)
endif()


# Build the netcdf c library.
set(NETCDF_LIBRARY "${PROJECT_BINARY_DIR}/lib/libnetcdf${LIB_SUFFIX}")

if (NOT EXISTS ${NETCDF_LIBRARY})
  message("Preparing netcdf-c...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/netcdf-c ${CMAKE_CURRENT_BINARY_DIR}/netcdf-c
                  RESULT_VARIABLE stat)

  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()
  message("Configuring netcdf-c library...")

#  SET(ZLIB_INCLUDE_DIRS "${PROJECT_BINARY_DIR}/include")
  message(STATUS "${CMAKE_COMMAND} ${ECOSIM_HDF5_LIBRARY} ${ECOSIM_HDF5_INCLUDE_DIR} ${ECOSIM_HDF5_HL_LIBRARY}")

  execute_process(COMMAND autoreconf -if WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf-c)

  set(NETCDF_C_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR})
  set(NETCDF_C_CONFIG_OPTS ${NETCDF_C_CONFIG_OPTS}  --disable-dap)
  set(NETCDF_C_CONFIG_OPTS ${NETCDF_C_CONFIG_OPTS}  --disable-parallel4)
  set(NETCDF_C_CONFIG_OPTS ${NETCDF_C_CONFIG_OPTS}  --enable-hdf5)
  if(BUILD_SHARED_LIBS)
    set(NETCDF_C_CONFIG_OPTS ${NETCDF_C_CONFIG_OPTS} --enable-shared=yes)
  else()
    set(NETCDF_C_CONFIG_OPTS ${NETCDF_C_CONFIG_OPTS} --enable-shared=no)
  endif()

  message("Netcdf-c opts: ${NETCDF_C_CONFIG_OPTS}")

  message("NETCDF_CPPFLAGS: ${NETCDF_CPPFLAGS}")
  message("NETCDF_LIBFLAGS: ${NETCDF_LIBFLAGS}")

  #Original version
  execute_process(COMMAND env CC=${CMAKE_C_COMPILER} env FC=${CMAKE_Fortran_COMPILER}
    env CXX=${CMAKE_CXX_COMPILER} env CXXFLAGS=${CMAKE_CXX_FLAGS} env CPPFLAGS=-I${PROJECT_BINARY_DIR}/include
    env CFLAGS=${CMAKE_C_FLAGS} env LDFLAGS=-L${PROJECT_BINARY_DIR}/lib
    ./configure ${NETCDF_C_CONFIG_OPTS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf-c
    OUTPUT_VARIABLE netcdf_c_config_log ERROR_VARIABLE netcdf_c_config_err
    RESULT_VARIABLE stat)

  #Attempt to change (might have made things worse)
  #execute_process(COMMAND env CC=${CMAKE_C_COMPILER} env FC=${CMAKE_Fortran_COMPILER}
  #  env CXX=${CMAKE_CXX_COMPILER} env CXXFLAGS=${CMAKE_CXX_FLAGS} env CPPFLAGS=${NETCDF_CPPFLAGS}
  #  env CFLAGS=${CMAKE_C_FLAGS} env LDFLAGS=${NETCDF_LIBFLAGS}
  #  env LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/lib/ env LIBS=${NETCDF_LIBFLAGS} 
  #  ./configure ${NETCDF_C_CONFIG_OPTS}
  #  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf-c
  #  OUTPUT_VARIABLE netcdf_c_config_log ERROR_VARIABLE netcdf_c_config_err
  #  RESULT_VARIABLE stat)

  message("In netcdf c build:")
  message("NETCDF_C_CONFIG_OPTS: ${NETCDF_C_CONFIG_OPTS}")

  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_c_config.log ${netcdf_c_config_log})
  if (NOT stat EQUAL 0)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_c_config_errors.log ${netcdf_c_config_err})

    if (TRAVIS_CI)
        message(WARNING "Configuring netcdf-c failed:")
        message(${netcdf_c_config_log})
        message(FATAL_ERROR ${netcdf_c_config_err})
      else()
        message(FATAL_ERROR "Error configuring netcdf-c library. See ${CMAKE_CURRENT_BINARY_DIR}/netcdf_c_config*.log for details.")
      endif()
  endif()

  message("Building and installing netcdf-c library...")
  execute_process(COMMAND make install -j1
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf-c
                  OUTPUT_VARIABLE netcdf_build_log ERROR_VARIABLE netcdf_build_err
                  RESULT_VARIABLE stat)

  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_c_build.log ${netcdf_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Error building netcdf-c:")
      message(${netcdf_c_build_log})
      message(FATAL_ERROR ${netcdf_c_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_c_build_errors.log ${netcdf_build_err})
      message(FATAL_ERROR "Error building netcdf-c library. See ${CMAKE_CURRENT_BINARY_DIR}/netcdf_c_build*.log for details.")
    endif()
  endif()

  # On Macs we need to change the rpath of shared libraries so we can install them elsewhere.
  if (APPLE AND BUILD_SHARED_LIBS)
    execute_process(COMMAND install_name_tool -id @rpath/libnetcdf${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/libnetcdf${LIB_SUFFIX})
  endif()
 #file(GLOB MY_PUBLIC_HEADERS "${PROJECT_BINARY_DIR}/OFF/netcdf*")
 #file(COPY ${MY_PUBLIC_HEADERS} DESTINATION ${PROJECT_BINARY_DIR}/include)
 #file(GLOB MY_PUBLIC_HEADERS "${PROJECT_BINARY_DIR}/OFF/type*")
 #file(COPY ${MY_PUBLIC_HEADERS} DESTINATION ${PROJECT_BINARY_DIR}/include)
endif()

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  if (${NETCDF_LIBRARY} STREQUAL "${PROJECT_BINARY_DIR}/lib/libnetcdf${LIB_SUFFIX}")
    file(GLOB nc_libs "${PROJECT_BINARY_DIR}/lib/libnetcdf*${LIB_SUFFIX}*")
    install(FILES ${nc_libs} DESTINATION lib)
    foreach(nc_util nc-config;nccopy;ncdump;ncgen;ncgen3)
      list(APPEND nc_utils ${PROJECT_BINARY_DIR}/bin/${nc_util})
    endforeach()
    install(FILES ${nc_utils} DESTINATION bin
            PERMISSIONS OWNER_READ;OWNER_WRITE;OWNER_EXECUTE;
                        GROUP_READ;GROUP_EXECUTE;WORLD_READ;WORLD_EXECUTE)
  
  message("NETCDF_LIBRARY: ${NETCDF_LIBRARY}")
  message("nc_libs: ${nc_libs}")
  message("nc_utils: ${nc_utils}")
  endif()
endif()

set(ECOSIM_TP_LIBS netcdf;${ECOSIM_TP_LIBS})

# Build the netcdf fortran library.
set(NETCDF_LIBRARY "${PROJECT_BINARY_DIR}/lib/libnetcdff${LIB_SUFFIX}")
if (NOT EXISTS ${NETCDF_LIBRARY})
  message("Preparing netcdf-fortran...")
  execute_process(COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/netcdf-fortran ${CMAKE_CURRENT_BINARY_DIR}/netcdf-fortran
                  RESULT_VARIABLE stat)
  if (NOT stat EQUAL 0)
    message(FATAL_ERROR "Copying library source failed.")
  endif()

  message("Configuring netcdf-fortran library...")

  execute_process(COMMAND autoreconf -if WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf-fortran)

  set(NETCDF_F_CONFIG_OPTS --prefix=${PROJECT_BINARY_DIR})
  if(BUILD_SHARED_LIBS)
    set(NETCDF_F_CONFIG_OPTS ${NETCDF_F_CONFIG_OPTS} --enable-shared=yes)
  else()
    set(NETCDF_F_CONFIG_OPTS ${NETCDF_F_CONFIG_OPTS} --enable-shared=no)
  endif()

  set(NETCDF_CPPFLAGS -I${PROJECT_BINARY_DIR}/include/)

  set(NETCDF_LIBFLAGS "")
  set(RESULT_VAR 0)
  execute_process(
    COMMAND ${PROJECT_BINARY_DIR}/bin/nc-config --libs
    OUTPUT_VARIABLE NETCDF_LIBFLAGS
    RESULT_VARIABLE RESULT_VAR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if (RESULT_VAR EQUAL 0)
    message("NETCDF_LIBFLAGS: ${NETCDF_LIBFLAGS}")
else ()
    message("Command failed with result: ${RESULT_VAR}")
endif ()

#  set(NETCDF_LIBFLAGS -L${PROJECT_BINARY_DIR}/lib/)

message("Netcdf-f opts: ${NETCDF_F_CONFIG_OPTS}")

  execute_process(COMMAND env CC=${CMAKE_C_COMPILER} env CFLAGS=${CMAKE_C_FLAGS}
     env CPPFLAGS=${NETCDF_CPPFLAGS} env LDFLAGS=${NETCDF_LIBFLAGS} env FC=${CMAKE_Fortran_COMPILER}
     env LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}/lib/ 
     env LIBS=${NETCDF_LIBFLAGS} ./configure  ${NETCDF_F_CONFIG_OPTS}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf-fortran
                  OUTPUT_VARIABLE netcdf_fortran_config_log ERROR_VARIABLE netcdf_fortran_config_err
                  RESULT_VARIABLE stat)

  message("In netcdf f build:")
  message("NETCDF_C_CONFIG_OPTS: ${NETCDF_F_CONFIG_OPTS}")

  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_fortran_config.log ${netcdf_fortran_config_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Error configuring netcdf-fortran:")
      message(${netcdf_fortran_config_log})
      message(FATAL_ERROR ${netcdf_fortran_config_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_fortran_config_errors.log ${netcdf_fortran_config_err})
      message(FATAL_ERROR "Error configuring netcdf-fortran library. See ${CMAKE_CURRENT_BINARY_DIR}/netcdf_fortran_config*.log for details.")
    endif()
  endif()

  message("Building and installing netcdf-fortran library...")
#  execute_process(COMMAND "${CMAKE_COMMAND}"
#    ${NETCDF_FORTRAN_CMAKE_OPTS}
#    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf-fortran
#    OUTPUT_VARIABLE netcdf_fortran_config_log ERROR_VARIABLE netcdf_fortran_config_err
#                   RESULT_VARIABLE stat)

  execute_process(COMMAND make install -j4
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/netcdf-fortran
                  OUTPUT_VARIABLE netcdf_fortran_build_log ERROR_VARIABLE netcdf_fortran_build_err
                  RESULT_VARIABLE stat)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_fortran_build.log ${netcdf_fortran_build_log})
  if (NOT stat EQUAL 0)
    if (TRAVIS_CI)
      message(WARNING "Error building netcdf-fortran")
      message(${netcdf_fortran_build_log})
      message(FATAL_ERROR ${netcdf_fortran_build_err})
    else()
      file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/netcdf_fortran_build_errors.log ${netcdf_fortran_build_err})
      message(FATAL_ERROR "Error building netcdf-fortran library. See ${CMAKE_CURRENT_BINARY_DIR}/netcdf_fortran_build*.log for details.")
    endif()
  endif()

  # On Macs we need to change the rpath of shared libraries so we can install them elsewhere.
  if (APPLE AND BUILD_SHARED_LIBS)
    execute_process(COMMAND install_name_tool -id @rpath/libnetcdff${LIB_SUFFIX} ${PROJECT_BINARY_DIR}/lib/libnetcdff${LIB_SUFFIX})
  endif()
 file(GLOB MY_PUBLIC_HEADERS "${PROJECT_BINARY_DIR}/OFF/netcdf*")
 file(COPY ${MY_PUBLIC_HEADERS} DESTINATION ${PROJECT_BINARY_DIR}/include)
 file(GLOB MY_PUBLIC_HEADERS "${PROJECT_BINARY_DIR}/OFF/type*")
 file(COPY ${MY_PUBLIC_HEADERS} DESTINATION ${PROJECT_BINARY_DIR}/include)
endif()

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  if (${NETCDF_LIBRARY} STREQUAL "${PROJECT_BINARY_DIR}/lib/libnetcdff${LIB_SUFFIX}")
    file(GLOB nc_libs "${PROJECT_BINARY_DIR}/lib/libnetcdff*${LIB_SUFFIX}*")
    install(FILES ${nc_libs} DESTINATION lib)
    foreach(nc_util nf-config)
      list(APPEND nc_utils ${PROJECT_BINARY_DIR}/bin/${nc_util})
    endforeach()
    install(FILES ${nc_utils} DESTINATION bin
            PERMISSIONS OWNER_READ;OWNER_WRITE;OWNER_EXECUTE;
                        GROUP_READ;GROUP_EXECUTE;WORLD_READ;WORLD_EXECUTE)
  message("NETCDF_LIBRARY: ${NETCDF_LIBRARY}")
  message("nc_libs: ${nc_libs}")
  message("nc_utils: ${nc_utils}")
  endif()
endif()
set(ECOSIM_TP_LIBS netcdff;${ECOSIM_TP_LIBS})

# Add all the libraries to the build system at large.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ECOSIM_TP_C_FLAGS}" PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ECOSIM_TP_Fortran_FLAGS}" PARENT_SCOPE)
#set(ECOSIM_LIBRARIES ${ECOSIM_LIBRARIES};${ECOSIM_TP_LIBS} PARENT_SCOPE)
#set(ECOSIM_INCLUDE_DIRS ${ECOSIM_INCLUDE_DIRS};${ECOSIM_TP_INCDIRS} PARENT_SCOPE)

set(NETCDF_LIBFLAGS "")
set(RESULT_VAR 0)
execute_process(
  COMMAND ${PROJECT_BINARY_DIR}/bin/nc-config --prefix
  OUTPUT_VARIABLE TPL_NETCDF_INCLUDE_DIRS
  RESULT_VARIABLE RESULT_VAR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(TPL_NETCDF_LIBRARIES "")
set(RESULT_VAR 0)
execute_process(
  COMMAND ${PROJECT_BINARY_DIR}/bin/nf-config --flibs
  OUTPUT_VARIABLE TPL_NETCDF_LIBRARIES
  RESULT_VARIABLE RESULT_VAR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)


set(TPL_NETCDF_LIBRARIES ${TPL_NETCDF_LIBRARIES} PARENT_SCOPE)
set(TPL_NETCDF_INCLUDE_DIRS "${TPL_NETCDF_INCLUDE_DIRS}/include" PARENT_SCOPE)
