function(add_ecosim_library lib)
  add_library(${lib} ${ARGN})
  if (BUILD_SHARED_LIBS)
	  target_link_libraries(${lib} ${ECOSIM_LIBRARIES})
  endif()
endfunction(add_ecosim_library)

