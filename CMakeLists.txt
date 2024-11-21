project(VoroPlusPlus_JL_wrapper)

cmake_minimum_required(VERSION 3.5)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

find_package(JlCxx)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")

find_package(VORO)
get_target_property(VORO_location VORO::voro++ LOCATION)
get_filename_component(VORO_location ${VORO_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location};${VORO_location}")

message(STATUS "Found Voro++ at ${VORO_location}")

add_library(voro++wrap SHARED src/wrapper.cpp)

target_link_libraries(voro++wrap JlCxx::cxxwrap_julia VORO::voro++)

install(TARGETS
  voro++wrap
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION lib)
