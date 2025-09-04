<<<<<<< HEAD
set(TECIO_SHARED TRUE)
set(TECIO_MPI FALSE)
=======
set(TECIO_SHARED FALSE)
<<<<<<< HEAD
set(TECIO_MPI TRUE)
>>>>>>> fceb2c53 (ReadSubzoneTecplot: clean PR (only module + CMake))
=======
set(TECIO_MPI FALSE)
>>>>>>> 69a4cd23 (chore: remove MPI from CmakeLists)
include(tecio.cmake)

set(HEADERS ReadSubzoneTecplot.h)

set(SOURCES ReadSubzoneTecplot.cpp)

#add_definitions(-DPYBIND11_DETAILED_ERROR_MESSAGES)
add_module(ReadSubzoneTecplot "read Sz Tecplot data " ${HEADERS} ${SOURCES})

<<<<<<< HEAD
<<<<<<< HEAD
target_link_libraries(ReadSubzoneTecplot PRIVATE vistle_tecio)
<<<<<<< HEAD
=======
target_link_libraries(ReadSubzoneTecplot vistle_tecio)
=======
target_link_libraries(ReadSubzoneTecplot PRIVATE vistle_tecio)
>>>>>>> 0740d5da (chore: apply pre-commit auto-fixes (EOF + clang-format), add missing scope to target libraries)

# If tecio.cmake does NOT already do this, add:
if(TECIO_MPI)
    find_package(MPI REQUIRED)
    target_compile_definitions(ReadSubzoneTecplot PRIVATE TECIOMPI)
    target_link_libraries(ReadSubzoneTecplot PRIVATE MPI::MPI_CXX)
endif()
>>>>>>> fceb2c53 (ReadSubzoneTecplot: clean PR (only module + CMake))
=======
>>>>>>> 69a4cd23 (chore: remove MPI from CmakeLists)
