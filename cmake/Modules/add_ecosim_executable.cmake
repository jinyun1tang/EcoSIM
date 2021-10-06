function(add_ecosim_executable exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${ECOSIM_LIBRARIES})
endfunction(add_ecosim_executable)

