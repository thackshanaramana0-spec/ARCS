cmake_minimum_required(VERSION 3.16)
project(arcs_cpp VERSION 2.2.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Compiler flags ────────────────────────────────────────────────────────────
option(ARCS_NATIVE    "Build with -march=native (fastest on this machine, not portable)" ON)
option(ARCS_PORTABLE  "Build with -msse4.1 (runs on any x86-64 with SSE4.1, ~2010+)" OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-O3 -funroll-loops -Wall -Wextra -Wno-unused-parameter)
    if(ARCS_PORTABLE)
        add_compile_options(-msse4.1)
    elseif(ARCS_NATIVE)
        add_compile_options(-march=native)
    endif()
elseif(MSVC)
    add_compile_options(/O2 /W4)
    if(ARCS_NATIVE)
        add_compile_options(/arch:AVX2)
    endif()
endif()

# ── OpenMP ────────────────────────────────────────────────────────────────────
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP ${OpenMP_CXX_VERSION} found — parallel KNN enabled")
else()
    message(STATUS "OpenMP NOT found — single-threaded KNN")
endif()

# ── Dependencies ──────────────────────────────────────────────────────────────
find_package(ZLIB REQUIRED)

# LZMA — try standard CMake module first (Linux/macOS), then MinGW paths (Windows)
find_package(LibLZMA QUIET)
if(LibLZMA_FOUND OR LIBLZMA_FOUND)
    message(STATUS "Found liblzma via find_package")
else()
    # Windows MinGW fallback
    find_library(LIBLZMA_LIB lzma HINTS C:/mingw64/lib)
    find_path(LIBLZMA_INC lzma.h HINTS C:/mingw64/include)
    if(LIBLZMA_LIB AND LIBLZMA_INC)
        set(LIBLZMA_FOUND TRUE)
        set(LIBLZMA_LIBRARIES ${LIBLZMA_LIB})
        set(LIBLZMA_INCLUDE_DIRS ${LIBLZMA_INC})
        message(STATUS "Found liblzma (MinGW): ${LIBLZMA_LIB}")
    else()
        message(WARNING "liblzma NOT found — LZMA disabled. Amplicon compression (DS-5) will fall back to zlib-9 and produce archives that do not match the paper's 1.224 MB result. Install: liblzma-dev (apt) or vcpkg install liblzma:x64-windows.")
    endif()
endif()

# ── zstd ──────────────────────────────────────────────────────────────────────
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(ZSTD QUIET libzstd)
endif()
if(NOT ZSTD_FOUND)
    find_library(ZSTD_LIB zstd
        HINTS
            "$ENV{HOME}/zstd_build/zstd-1.5.6/lib"
            C:/mingw64/lib /usr/lib/x86_64-linux-gnu /usr/lib /usr/local/lib)
    find_path(ZSTD_INC zstd.h
        HINTS
            "$ENV{HOME}/zstd_build/zstd-1.5.6/lib"
            C:/mingw64/include /usr/include /usr/local/include)
    if(ZSTD_LIB AND ZSTD_INC)
        set(ZSTD_FOUND TRUE)
        set(ZSTD_LIBRARIES ${ZSTD_LIB})
        set(ZSTD_INCLUDE_DIRS ${ZSTD_INC})
    endif()
endif()
if(ZSTD_FOUND)
    message(STATUS "Found zstd: ${ZSTD_LIBRARIES} — large blobs use zstd-19 (lower RAM)")
else()
    message(STATUS "zstd NOT found — large blobs use LZMA-9 (install libzstd-dev)")
endif()

# Helper to apply zstd to a target
function(target_link_zstd tgt)
    if(ZSTD_FOUND)
        target_compile_definitions(${tgt} PRIVATE ARCS_HAVE_ZSTD=1)
        if(ZSTD_INCLUDE_DIRS)
            target_include_directories(${tgt} PRIVATE ${ZSTD_INCLUDE_DIRS})
        endif()
        target_link_libraries(${tgt} PRIVATE ${ZSTD_LIBRARIES})
    endif()
endfunction()

# ── Helper to enable libsais's OpenMP-parallel SA/PLCP/LCP construction ───────
# (used by src/sa_apsp.cpp) — every target that compiles ARCS_SOURCES needs
# this applied consistently, since libsais.h itself guards the *_omp function
# declarations behind LIBSAIS_OPENMP: a target missing this define would fail
# to link against sa_apsp.cpp's calls to those functions.
function(target_enable_libsais_openmp tgt)
    if(OpenMP_C_FOUND AND OpenMP_CXX_FOUND)
        target_compile_definitions(${tgt} PRIVATE LIBSAIS_OPENMP=1)
        target_link_libraries(${tgt} PRIVATE OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
    else()
        message(WARNING "OpenMP not found for both C and C++ — libsais will run single-threaded for target ${tgt}")
    endif()
endfunction()

# ── libbsc (bundled) ──────────────────────────────────────────────────────────
set(LIBBSC_SOURCES
    third_party/libbsc/adler32/adler32.cpp
    third_party/libbsc/bwt/bwt.cpp
    third_party/libbsc/bwt/divsufsort/divsufsort.c
    third_party/libbsc/coder/coder.cpp
    third_party/libbsc/coder/qlfc/qlfc.cpp
    third_party/libbsc/coder/qlfc/qlfc_model.cpp
    third_party/libbsc/filters/detectors.cpp
    third_party/libbsc/filters/preprocessing.cpp
    third_party/libbsc/libbsc/libbsc.cpp
    third_party/libbsc/lzp/lzp.cpp
    third_party/libbsc/platform/platform.cpp
)

# ── Sources ───────────────────────────────────────────────────────────────────
set(ARCS_SOURCES
    src/fastq_io.cpp
    src/kmer.cpp
    src/debruijn.cpp
    src/minimizer.cpp
    src/mapper.cpp
    src/hamming_mst.cpp
    src/mst_encoder.cpp
    src/rans_model.cpp
    src/position_enc.cpp
    src/akc.cpp
    src/container.cpp
    src/chain_encoder.cpp
    src/sa_apsp.cpp
    src/hash_apsp.cpp
    src/fm_apsp.cpp
    src/vodbg_pg.cpp
    src/kmer_anchor_pg.cpp
    src/dna_coder.cpp
    src/repeat_elim.cpp
    src/qual_cm.cpp
    src/encoder.cpp
    src/decoder.cpp
    src/caller.cpp
    src/tensor.cpp
    src/bsc_codec.cpp
    ${LIBBSC_SOURCES}
    third_party/libsais/libsais.c
    third_party/libsais/libsais64.c
)

# ── Main binary ───────────────────────────────────────────────────────────────
add_executable(arcs src/main.cpp ${ARCS_SOURCES})

target_include_directories(arcs PRIVATE
    src/
    third_party/
    third_party/libbsc
)

find_package(Threads REQUIRED)
target_link_libraries(arcs PRIVATE ZLIB::ZLIB Threads::Threads)
if(OpenMP_CXX_FOUND)
    target_link_libraries(arcs PRIVATE OpenMP::OpenMP_CXX)
endif()
target_enable_libsais_openmp(arcs)

if(LIBLZMA_FOUND)
    target_compile_definitions(arcs PRIVATE ARCS_HAVE_LZMA=1)
    target_include_directories(arcs PRIVATE ${LIBLZMA_INCLUDE_DIRS})
    target_link_libraries(arcs PRIVATE ${LIBLZMA_LIBRARIES})
endif()
target_link_zstd(arcs)

# ── Helper to apply LZMA to a target ──────────────────────────────────────────
function(target_link_lzma tgt)
    if(LIBLZMA_FOUND)
        target_compile_definitions(${tgt} PRIVATE ARCS_HAVE_LZMA=1)
        target_include_directories(${tgt} PRIVATE ${LIBLZMA_INCLUDE_DIRS})
        target_link_libraries(${tgt} PRIVATE ${LIBLZMA_LIBRARIES})
    endif()
endfunction()


# ── Tests ─────────────────────────────────────────────────────────────────────
enable_testing()

add_executable(test_kmer tests/test_kmer.cpp ${ARCS_SOURCES})
target_include_directories(test_kmer PRIVATE src/ third_party/)
target_link_libraries(test_kmer PRIVATE ZLIB::ZLIB)
target_link_lzma(test_kmer)
target_link_zstd(test_kmer)
target_enable_libsais_openmp(test_kmer)
add_test(NAME kmer COMMAND test_kmer)

add_executable(test_debruijn tests/test_debruijn.cpp ${ARCS_SOURCES})
target_include_directories(test_debruijn PRIVATE src/ third_party/)
target_link_libraries(test_debruijn PRIVATE ZLIB::ZLIB)
target_link_lzma(test_debruijn)
target_link_zstd(test_debruijn)
target_enable_libsais_openmp(test_debruijn)
add_test(NAME debruijn COMMAND test_debruijn)

add_executable(test_rans tests/test_rans.cpp ${ARCS_SOURCES})
target_include_directories(test_rans PRIVATE src/ third_party/)
target_link_libraries(test_rans PRIVATE ZLIB::ZLIB)
target_link_lzma(test_rans)
target_link_zstd(test_rans)
target_enable_libsais_openmp(test_rans)
add_test(NAME rans COMMAND test_rans)

add_executable(test_mapper tests/test_mapper.cpp ${ARCS_SOURCES})
target_include_directories(test_mapper PRIVATE src/ third_party/)
target_link_libraries(test_mapper PRIVATE ZLIB::ZLIB)
target_link_lzma(test_mapper)
target_link_zstd(test_mapper)
target_enable_libsais_openmp(test_mapper)
add_test(NAME mapper COMMAND test_mapper)

add_executable(test_roundtrip tests/test_roundtrip.cpp ${ARCS_SOURCES})
target_include_directories(test_roundtrip PRIVATE src/ third_party/)
target_link_libraries(test_roundtrip PRIVATE ZLIB::ZLIB)
target_link_lzma(test_roundtrip)
target_link_zstd(test_roundtrip)
target_enable_libsais_openmp(test_roundtrip)
add_test(NAME roundtrip COMMAND test_roundtrip)

add_executable(test_qual_rans tests/test_qual_rans.cpp ${ARCS_SOURCES})
target_include_directories(test_qual_rans PRIVATE src/ third_party/)
target_link_libraries(test_qual_rans PRIVATE ZLIB::ZLIB)
target_link_lzma(test_qual_rans)
target_link_zstd(test_qual_rans)
target_enable_libsais_openmp(test_qual_rans)
add_test(NAME qual_rans COMMAND test_qual_rans)

add_executable(test_dna_coder tests/test_dna_coder.cpp ${ARCS_SOURCES})
target_include_directories(test_dna_coder PRIVATE src/ third_party/)
target_link_libraries(test_dna_coder PRIVATE ZLIB::ZLIB)
target_link_lzma(test_dna_coder)
target_link_zstd(test_dna_coder)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_dna_coder PRIVATE OpenMP::OpenMP_CXX)
endif()
target_enable_libsais_openmp(test_dna_coder)
add_test(NAME dna_coder COMMAND test_dna_coder)

add_executable(test_vodbg_calldata tests/test_vodbg_calldata.cpp ${ARCS_SOURCES})
target_include_directories(test_vodbg_calldata PRIVATE src/ third_party/)
target_link_libraries(test_vodbg_calldata PRIVATE ZLIB::ZLIB)
target_link_lzma(test_vodbg_calldata)
target_link_zstd(test_vodbg_calldata)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_vodbg_calldata PRIVATE OpenMP::OpenMP_CXX)
endif()
target_enable_libsais_openmp(test_vodbg_calldata)
add_test(NAME vodbg_calldata COMMAND test_vodbg_calldata)

# test_sa_large_scale: NOT wired into `add_test` deliberately — it builds a
# >2.147B-char (~2.3GB) text and needs ~20GB peak RAM + real wall time to
# construct+verify, which would slow down the fast default ctest suite (and
# violate the "7/7 ctests before any benchmark phase" quick-check budget).
# Run explicitly (./build/test_sa_large_scale) when validating the
# >INT32_MAX/libsais64 code path specifically, not as part of routine ctest.
add_executable(test_sa_large_scale tests/test_sa_large_scale.cpp src/sa_apsp.cpp third_party/libsais/libsais.c third_party/libsais/libsais64.c)
target_include_directories(test_sa_large_scale PRIVATE src/ third_party/)
target_enable_libsais_openmp(test_sa_large_scale)

# ── Install ───────────────────────────────────────────────────────────────────
install(TARGETS arcs RUNTIME DESTINATION bin)
