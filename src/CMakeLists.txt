set(MODULE_NAME compMath)

message("comp_math_build " ${MODULE_NAME})

file(GLOB_RECURSE SOURCES *.h *.cpp *.hpp)

add_library(${MODULE_NAME} STATIC ${SOURCES})

set_target_properties(${MODULE_NAME} PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(${MODULE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(${MODULE_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

target_link_libraries(${MODULE_NAME} PUBLIC eigen-3.4.0)
