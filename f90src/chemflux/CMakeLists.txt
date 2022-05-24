add_subdirectory(nonsalt)
add_subdirectory(salt)
set(ECOSIM_LIBRARIES nonsalt;salt;${ECOSIM_LIBRARIES} PARENT_SCOPE)
