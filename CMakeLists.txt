cmake_minimum_required(VERSION 3.0)
project(LIV_neutrino)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2 ")
include_directories(/usr/local/include)
# mpi test 
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})




set(target_sources reactor_fluxcheck.c myio.c Lorentz_violation_lat.h)
list(TRANSFORM target_sources PREPEND "code_source/")
add_executable(reactor_fluxcheck ${target_sources})

target_link_libraries(reactor_fluxcheck
    -L/usr/local/lib -Wl,-rpath -Wl,/usr/local/lib -lglobes -L/usr/local/lib
    -lgsl -lgslcblas -lm -lm -static
)

set(target_sources LIV_minizer_lat.c myio.c DoubleChooz_chi.h RENO_chi.h DayaBay_chi.h Lorentz_violation_lat.h glbmath_lat.h)
list(TRANSFORM target_sources PREPEND "code_source/")
add_executable(LIV_minizer_lat ${target_sources})

target_link_libraries(LIV_minizer_lat
    -L/usr/local/lib -Wl,-rpath -Wl,/usr/local/lib -lglobes -L/usr/local/lib
    -lgsl -lgslcblas -lm -lm -static
)