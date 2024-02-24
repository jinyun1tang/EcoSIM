function(add_ecosim_library lib)
  
  if (BUILD_SHARED_LIBS)
    add_library(${lib} ${ARGN})
  else()
    add_library(${lib} STATIC ${ARGN})    
  endif()
  target_link_libraries(${lib} ${ECOSIM_LIBRARIES} ${ECOSIM_TPLS})  
endfunction(add_ecosim_library)

