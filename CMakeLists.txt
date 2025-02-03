cmake_minimum_required(VERSION 3.0)
project(LIV_neutrino)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2 ")
include_directories(/usr/local/include)
# mpi test 
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})




set(target_sources reactor_fluxcheck.c myio.c)
list(TRANSFORM target_sources PREPEND "code_source/")
add_executable(reactor_fluxcheck ${target_sources})

target_link_libraries(reactor_fluxcheck
    -L/usr/local/lib -Wl,-rpath -Wl,/usr/local/lib -lglobes -L/usr/local/lib
    -lgsl -lgslcblas -lm -lm -static
)
