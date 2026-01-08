cmake_minimum_required(VERSION 3.16)

project(MESSI_SFA C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(MESSI_ENABLE_OPENMP "Enable OpenMP if available." ON)
option(MESSI_ENABLE_SIMD "Enable SIMD flags if supported by the compiler." ON)
option(MESSI_ENABLE_NATIVE_ARCH "Enable -march=native for local builds." ON)
option(MESSI_REQUIRE_FFTW "Fail configure if FFTW3f is missing." ON)
option(MESSI_REQUIRE_READLINE "Fail configure if readline is missing." ON)

set(ADS_SOURCES
  src/ads/isax_file_loaders.c
  src/ads/isax_first_buffer_layer.c
  src/ads/isax_index.c
  src/ads/isax_node.c
  src/ads/isax_node_buffer.c
  src/ads/isax_node_record.c
  src/ads/isax_node_split.c
  src/ads/isax_query_engine.c
  src/ads/isax_visualize_index.c
  src/ads/pqueue.c
  src/ads/sax/sax.c
  src/ads/sax/ts.c
  src/ads/inmemory_query_engine.c
  src/ads/parallel_inmemory_query_engine.c
  src/ads/inmemory_index_engine.c
  src/ads/parallel_index_engine.c
  src/ads/parallel_query_engine.c
  src/ads/inmemory_topk_engine.c
  src/ads/sfa/dft.c
  src/ads/sfa/sfa.c
  src/ads/spartan/pca.c
  src/ads/spartan/spartan.c
  src/ads/calc_utils.c
)

add_library(ads STATIC ${ADS_SOURCES})
target_include_directories(ads
  PUBLIC
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/include
)
target_compile_options(ads PRIVATE -fcommon)

if(MESSI_ENABLE_NATIVE_ARCH)
  target_compile_options(ads PRIVATE -march=native)
endif()

if(MESSI_ENABLE_SIMD)
  include(CheckCCompilerFlag)
  check_c_compiler_flag("-mavx2" MESSI_HAVE_AVX2)
  if(MESSI_HAVE_AVX2)
    target_compile_options(ads PRIVATE -mavx -mavx2 -msse3)
  endif()
endif()

if(MESSI_ENABLE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_C_FOUND)
    target_link_libraries(ads PRIVATE OpenMP::OpenMP_C)
  else()
    message(WARNING "OpenMP not found; building without -fopenmp.")
  endif()
endif()

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(FFTW3F QUIET IMPORTED_TARGET fftw3f)
endif()

if(TARGET PkgConfig::FFTW3F)
  set(MESSI_FFTW_TARGET PkgConfig::FFTW3F)
else()
  set(FFTW_ROOT_HINTS
    ${FFTW_ROOT}
    $ENV{FFTW_ROOT}
  )
  find_path(FFTW3F_INCLUDE_DIR fftw3.h
    HINTS ${FFTW_ROOT_HINTS}
    PATH_SUFFIXES include
  )
  find_library(FFTW3F_LIBRARY fftw3f
    HINTS ${FFTW_ROOT_HINTS}
    PATH_SUFFIXES lib
  )

  if(FFTW3F_INCLUDE_DIR AND FFTW3F_LIBRARY)
    add_library(FFTW3F::FFTW3F UNKNOWN IMPORTED)
    set_target_properties(FFTW3F::FFTW3F PROPERTIES
      IMPORTED_LOCATION "${FFTW3F_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW3F_INCLUDE_DIR}"
    )
    set(MESSI_FFTW_TARGET FFTW3F::FFTW3F)
  endif()
endif()

if(NOT MESSI_FFTW_TARGET)
  if(MESSI_REQUIRE_FFTW)
    message(FATAL_ERROR "FFTW3f not found. Install fftw3f or set FFTW_ROOT.")
  else()
    message(WARNING "FFTW3f not found; build will likely fail without it.")
  endif()
endif()

add_executable(MESSI src/utils/MESSI.c)
target_include_directories(MESSI
  PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/include
)
target_compile_options(MESSI PRIVATE -fcommon)

if(MESSI_ENABLE_NATIVE_ARCH)
  target_compile_options(MESSI PRIVATE -march=native)
endif()

if(MESSI_ENABLE_SIMD AND MESSI_HAVE_AVX2)
  target_compile_options(MESSI PRIVATE -mavx -mavx2 -msse3)
endif()

if(MESSI_ENABLE_OPENMP AND OpenMP_C_FOUND)
  target_link_libraries(MESSI PRIVATE OpenMP::OpenMP_C)
endif()

target_link_libraries(MESSI PRIVATE ads)

if(MESSI_FFTW_TARGET)
  target_link_libraries(MESSI PRIVATE ${MESSI_FFTW_TARGET})
endif()

find_library(READLINE_LIBRARY readline)
if(READLINE_LIBRARY)
  target_link_libraries(MESSI PRIVATE ${READLINE_LIBRARY})
elseif(MESSI_REQUIRE_READLINE)
  message(FATAL_ERROR "readline not found. Install readline or set MESSI_REQUIRE_READLINE=OFF.")
else()
  message(WARNING "readline not found; build will likely fail without it.")
endif()

find_package(Threads REQUIRED)
target_link_libraries(MESSI PRIVATE Threads::Threads)

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
  target_link_libraries(MESSI PRIVATE ${MATH_LIBRARY})
endif()
