add_subdirectory(box_model)
add_subdirectory(layers_model)
set(ECOSIM_LIBRARIES box_model;layers_model;${ECOSIM_LIBRARIES} PARENT_SCOPE)
