function(add_ecosim_executable exe)
  add_executable(${exe} ${ARGN})

  message("We are building exe ${exe}, with:")
  message(STATUS "ECOSIM_LIBRARIES: ${ECOSIM_LIBRARIES}")
  message(STATUS "ECOSIM_TPLS: ${ECOSIM_TPLS}")

  target_link_libraries(${exe} ${ECOSIM_LIBRARIES} ${ECOSIM_TPLS})
endfunction(add_ecosim_executable)

