########################################
# Set up how to compile the source files
########################################

# Define the executable name
SET(2DPMEEXE 2dpme.exe)

# Add the source files
SET(2DPME_Sources ${SRC}/main.f90)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${2DPMEEXE} ${2DPME_Sources})

# This links foo to the bar library
#TARGET_LINK_LIBRARIES(${2DPMEEXE} ${SOME_LIBRARY})

#####################################
# Tell how to install this executable
#####################################

set_target_properties( ${2DPMEEXE}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
