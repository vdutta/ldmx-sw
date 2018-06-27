# declare playtest module
module(
  NAME playtest  
  DEPENDENCIES Framework Event DetDescr
  EXTERNAL_DEPENDENCIES ROOT Python
)

# generate ROOT dictionary for Event classes
root_generate_dictionary(HcalTrackDict ${Event_INCLUDE_DIR}/playtest/HcalTrack.h
MODULE ${PROJECT_NAME} LINKDEF ${Event_INCLUDE_DIR}/playtest/HcalTrackLink.h)

# install ROOT pcm file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${PROJECT_NAME}_rdict.pcm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
