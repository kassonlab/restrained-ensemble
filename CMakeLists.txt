cmake_minimum_required(VERSION 3.7)
project(Roux)

set(CMAKE_CXX_STANDARD 11)

#set(BOOST_ROOT /opt/boost/intel/openmpi_ib/)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.5 COMPONENTS mpi filesystem program_options regex REQUIRED)

if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIR})
endif ()

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(SOURCE_FILES main.cpp ensemble.cpp rouxfileio.cpp rouxalgorithms.cpp jsoncpp.cpp logging.cpp logging.h)
add_executable(Roux ${SOURCE_FILES})
if (MPI_COMPILE_FLAGS)
    set_target_properties(Roux PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set_target_properties(Roux PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif ()
include_directories(/home/jmh/gromacs-roux/src/programs/mdrun)
#include_directories(/home/jmhays/gromacs-roux-fb/src/programs/mdrun)
target_link_libraries(Roux ${Boost_LIBRARIES} ${MPI_LIBRARIES}
        /home/jmh/gromacs-roux/mpibuild/lib/libmdrun_lib.so
        /home/jmh/gromacs-roux/mpibuild/lib/libgromacs_mpi.so)
#        /home/jmhays/gromacs-roux-fb/mpibuild/lib/libmdrun_lib.so
#        /home/jmhays/gromacs-roux-fb/mpibuild/lib/libgromacs_mpi.so)

