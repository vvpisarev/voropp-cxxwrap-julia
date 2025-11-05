# Using the given filepath, load a version.json file into VERSION_MAJOR, VERSION_MINOR, and VERSION_PATCH cache variables
# As a macro this will be executed in the context of the CMakeLists that calls it.
macro(load_version_json version_filepath)
    if(EXISTS ${version_filepath})
        file(READ ${version_filepath} version_contents)
        # Note: variables are set to "keyname-NOTFOUND" if they are not present in the JSON string.
        # An error is emitted when the key is not found. This error is again emitted by the `project`
        # directive if you attempt to use a string like "<key>-NOTFOUND" as a version number. This
        # error will halt the build
        string(JSON VERSION_MAJOR GET ${version_contents} "major")
        string(JSON VERSION_MINOR GET ${version_contents} "minor")
        string(JSON VERSION_PATCH GET ${version_contents} "patch")
        set(SEMVER_STRING "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}" CACHE STRING "Semantic Version String for the project")
    else()
        message(FATAL_ERROR "Unable to load version info file at ${version_filepath}")
    endif()
endmacro()

if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

load_version_json(${CMAKE_CURRENT_LIST_DIR}/version.json)
project(voropp-cxxwrap-julia
    VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
    LANGUAGES C CXX
)

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

target_link_libraries(voro++wrap PUBLIC JlCxx::cxxwrap_julia JlCxx::cxxwrap_julia_stl VORO::voro++)

install(TARGETS
  voro++wrap
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION lib)
