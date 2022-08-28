add_subdirectory(box_chem)
add_subdirectory(layers_chem)
set(ECOSIM_LIBRARIES box_chem;layers_chem;${ECOSIM_LIBRARIES} PARENT_SCOPE)
