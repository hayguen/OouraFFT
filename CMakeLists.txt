cmake_minimum_required(VERSION 2.8)
project(OOURA_FFT)

option(OOURA_USE_FAST_MATH  "Build with fast math - non IEEE compliant" ON)
option(OOURA_USE_FLOAT_PREC "Build all libraries for single precision float - else: double" OFF)

if (WIN32 AND MINGW)
  option(OOURA_STATIC_BUILD "Build full static libraries on MinGW/Win32" ON)
else()
  set(OOURA_STATIC_BUILD false)
endif()

option(DISABLE_LINK_WITH_M "Disables linking with m library to build with clangCL from MSVC" OFF)

# C90 requires the gcc extensions for function attributes like always_inline
# C99 provides the function attributes: no gcc extensions required
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)

# detect compiler
set(COMPILER_IS_GCC false)
set(COMPILER_IS_CLANG false)
set(COMPILER_IS_MSVC false)
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  # using GCC
  set(COMPILER_IS_GCC true)
  message(STATUS "detected compiler: GCC")
  set(ARCH "native" CACHE STRING "gcc target architecture (-march): native/SSE:core2/AVX:sandybridge/ARM-NEON:armv7-a")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
  set(COMPILER_IS_CLANG true)
  message(STATUS "detected compiler: CLANG")
  set(ARCH "native" CACHE STRING "gcc target architecture (-march): native/SSE:core2/AVX:sandybridge")
elseif (CMAKE_C_COMPILER_ID STREQUAL "MSVC")  # others: "Clang", "Intel"
  set(COMPILER_IS_MSVC true)
  message(STATUS "detected compiler: MSVC")
  set(ARCH "AVX" CACHE STRING "msvc target architecture (/arch): SSE2/AVX")
else()
  message(WARNING "detected compiler: unknown!")
endif()


if(OOURA_STATIC_BUILD)
  # Special MINGW stuff here
  # see https://cmake.org/pipermail/cmake/2012-September/051970.html
  # see http://stackoverflow.com/questions/13768515/how-to-do-static-linking-of-libwinpthread-1-dll-in-mingw
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-libgcc -static")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++")
  set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_C_FLAGS} -static-libgcc -s")
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS} -static-libgcc -static-libstdc++ -s")
endif()



########################################################################
# select the release build type by default to get optimization flags
########################################################################
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
   message(STATUS "Build type not specified: defaulting to release.")
endif(NOT CMAKE_BUILD_TYPE)


# link with math library?
if(DISABLE_LINK_WITH_M)
  set(MATHLIB "")
else()
  message(STATUS "INFO: detected NO MSVC: ${CMAKE_C_COMPILER_ID}: will link math lib m")
  if (COMPILER_IS_MSVC)
    message(WARNING "INFO: detected MSVC: will not link math lib m")
    set(MATHLIB "")
  else()
    set(MATHLIB "m")
  endif()
endif()


# determine compiler specific options for target architecture
# and fast math
set(OPT_ARCH "")
set(OPT_FAST_MATH "")

if (COMPILER_IS_MSVC)
  if (ARCH STREQUAL "")
    message(WARNING "ARCH is not set!")
  else()
    message(STATUS "ARCH is ${ARCH}")
    set(OPT_ARCH "/arch:${ARCH}")
  endif()
  if (OOURA_USE_FAST_MATH)
    set(OPT_FAST_MATH "/fp:fast")
  endif()
endif()

if (COMPILER_IS_GCC OR COMPILER_IS_CLANG)
  if (ARCH STREQUAL "")
    message(WARNING "ARCH is not set!")
  else()
    message(STATUS "ARCH is ${ARCH}")
    set(OPT_ARCH "-march=${ARCH}")
  endif()
  if (OOURA_USE_FAST_MATH)
    set(OPT_FAST_MATH "-ffast-math")
  endif()
endif()

message(STATUS "compiler option for target architecture: ${OPT_ARCH}")
message(STATUS "compiler option for fast math: ${OPT_FAST_MATH}")

######################################################


# Simple Version I   (radix 4,2)
add_library(Ooura1dFFTradix4_h STATIC fft4g_h.c include/ooura/fftxg_h.h)

add_executable(test_radix4_h sample1/testxg_h.c)
target_link_libraries(test_radix4_h Ooura1dFFTradix4_h)

add_executable(bench_pi_fft_radix4_h sample2/pi_fft_h.c)
target_link_libraries(bench_pi_fft_radix4_h Ooura1dFFTradix4_h)

######################################################

# Fast Version I   (radix 4,2)
add_library(Ooura1dFFTradix4 STATIC fft4g.c include/ooura/fftxg.h)

add_executable(test_radix4 sample1/testxg.c)
target_link_libraries(test_radix4 Ooura1dFFTradix4)

add_executable(bench_pi_fft_radix4 sample2/pi_fft.c)
target_link_libraries(bench_pi_fft_radix4 Ooura1dFFTradix4)

######################################################

# Simple Version II  (radix 8,4,2)
add_library(Ooura1dFFTradix8_h STATIC fft8g_h.c include/ooura/fftxg_h.h)

add_executable(test_radix8_h sample1/testxg_h.c)
target_link_libraries(test_radix8_h Ooura1dFFTradix8_h)

add_executable(bench_pi_fft_radix8_h sample2/pi_fft_h.c)
target_link_libraries(bench_pi_fft_radix8_h Ooura1dFFTradix8_h)

######################################################

# Fast Version II  (radix 8,4,2)
add_library(Ooura1dFFTradix8 STATIC fft8g.c include/ooura/fftxg.h)

add_executable(test_radix8 sample1/testxg.c)
target_link_libraries(test_radix8 Ooura1dFFTradix8)

add_executable(bench_pi_fft_radix8 sample2/pi_fft.c)
target_link_libraries(bench_pi_fft_radix8 Ooura1dFFTradix8)

######################################################

# Simple Version III (Split-Radix) - no threading
add_library(Ooura1dFFTsplitradix_nt_h STATIC fftsg_h.c include/ooura/fftxg_h.h)

add_executable(test_splitradix_nt_h sample1/testxg_h.c)
target_link_libraries(test_splitradix_nt_h Ooura1dFFTsplitradix_nt_h)

add_executable(bench_pi_fft_splitradix_nt_h sample2/pi_fft_h.c)
target_link_libraries(bench_pi_fft_splitradix_nt_h Ooura1dFFTsplitradix_nt_h)

######################################################

# Fast Version III (Split-Radix) - no threading
add_library(Ooura1dFFTsplitradix_nt STATIC fftsg.c include/ooura/fftxg.h)

add_executable(test_splitradix_nt sample1/testxg.c)
target_link_libraries(test_splitradix_nt Ooura1dFFTsplitradix_nt)

add_executable(bench_pi_fft_splitradix_nt sample2/pi_fft.c)
target_link_libraries(bench_pi_fft_splitradix_nt Ooura1dFFTsplitradix_nt)

######################################################

# Simple Version III (Split-Radix) - multi threading (default): winthreads: 32K / 512K  else: pthreads: 8K / 64K
add_library(Ooura1dFFTsplitradix_mt_h STATIC fftsg_h.c include/ooura/fftxg_h.h)
if (WIN32)
  target_compile_definitions(Ooura1dFFTsplitradix_mt_h PRIVATE USE_CDFT_WINTHREADS)
else()
  target_compile_definitions(Ooura1dFFTsplitradix_mt_h PRIVATE USE_CDFT_PTHREADS)
  target_link_libraries(Ooura1dFFTsplitradix_mt_h pthread)
endif()

add_executable(test_splitradix_mt_h sample1/testxg_h.c)
target_link_libraries(test_splitradix_mt_h Ooura1dFFTsplitradix_mt_h)

add_executable(bench_pi_fft_splitradix_mt_h sample2/pi_fft_h.c)
target_link_libraries(bench_pi_fft_splitradix_mt_h Ooura1dFFTsplitradix_mt_h)

######################################################

# Fast Version III (Split-Radix) - multi threading (default): winthreads: 32K / 512K  else: pthreads: 8K / 64K
add_library(Ooura1dFFTsplitradix_mt STATIC fftsg.c include/ooura/fftxg.h)
if (WIN32)
  target_compile_definitions(Ooura1dFFTsplitradix_mt PRIVATE USE_CDFT_WINTHREADS)
else()
  target_compile_definitions(Ooura1dFFTsplitradix_mt PRIVATE USE_CDFT_PTHREADS)
  target_link_libraries(Ooura1dFFTsplitradix_mt pthread)
endif()

add_executable(test_splitradix_mt sample1/testxg.c)
target_link_libraries(test_splitradix_mt Ooura1dFFTsplitradix_mt)

add_executable(bench_pi_fft_splitradix_mt sample2/pi_fft.c)
target_link_libraries(bench_pi_fft_splitradix_mt Ooura1dFFTsplitradix_mt)

######################################################

set(OOURA_LIBS
  Ooura1dFFTradix4_h Ooura1dFFTradix4
  Ooura1dFFTradix8_h Ooura1dFFTradix8
  Ooura1dFFTsplitradix_nt_h Ooura1dFFTsplitradix_nt
  Ooura1dFFTsplitradix_mt_h Ooura1dFFTsplitradix_mt
)

set(OOURA_TESTS
  test_radix4_h test_radix4
  test_radix8_h test_radix8
  test_splitradix_nt_h test_splitradix_nt
  test_splitradix_mt_h test_splitradix_mt
)

set(OOURA_BENCHES
  bench_pi_fft_radix4_h bench_pi_fft_radix4
  bench_pi_fft_radix8_h bench_pi_fft_radix8
  bench_pi_fft_splitradix_nt_h bench_pi_fft_splitradix_nt
  bench_pi_fft_splitradix_mt_h bench_pi_fft_splitradix_mt
)


foreach(X IN LISTS OOURA_LIBS)
  # message(STATUS "set lib specs for ${X}")
  if (OOURA_USE_FLOAT_PREC)
    target_compile_definitions(${X} PUBLIC OOURA_SINGLE_PREC)
  endif()
  target_include_directories(${X} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
endforeach()

foreach(X IN LISTS OOURA_LIBS OOURA_TESTS OOURA_BENCHES)
  # message(STATUS "set specs for ${X}")
  target_compile_options(${X} PRIVATE $<$<C_COMPILER_ID:GNU>:-Wall -Wextra -pedantic>)
  target_compile_options(${X} PRIVATE "${OPT_ARCH}" "${OPT_FAST_MATH}")
  target_link_libraries(${X} ${MATHLIB})
endforeach()


######################################################

enable_testing()

add_test(NAME test_radix4_n512            COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_radix4" "512")
add_test(NAME test_radix8_n512            COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_radix8" "512")
add_test(NAME test_splitradix_nt_n512     COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_splitradix_nt" "512")
add_test(NAME test_splitradix_mt_n512     COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_splitradix_mt" "512")

add_test(NAME test_radix4_n65536          COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_radix4" "65536")
add_test(NAME test_radix8_n65536          COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_radix8" "65536")
add_test(NAME test_splitradix_nt_n65536   COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_splitradix_nt" "65536")
add_test(NAME test_splitradix_mt_n65536   COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_splitradix_mt" "65536")


add_test(NAME test_radix4_h_n512          COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_radix4_h" "512")
add_test(NAME test_radix8_h_n512          COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_radix8_h" "512")
add_test(NAME test_splitradix_nt_h_n512   COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_splitradix_nt_h" "512")
add_test(NAME test_splitradix_mt_h_n512   COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_splitradix_mt_h" "512")

add_test(NAME test_radix4_h_n65536        COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_radix4_h" "65536")
add_test(NAME test_radix8_h_n65536        COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_radix8_h" "65536")
add_test(NAME test_splitradix_nt_h_n65536 COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_splitradix_nt_h" "65536")
add_test(NAME test_splitradix_mt_h_n65536 COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_splitradix_mt_h" "65536")
