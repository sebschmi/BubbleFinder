cmake_minimum_required(VERSION 3.20)
project(sbSPQRnew LANGUAGES CXX)

# ──────────────────────────────────────────────────────────────────
# 0. Global compile flags – perf friendly
#    *  -g                       : full DWARF debug info
#    *  -fno-omit-frame-pointer : keep frame pointers for perf
# ──────────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(
    -g
    -fno-omit-frame-pointer
    # -fno-inline
)

# Release = fast but still symbolised
# set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2")

set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG    "-O2 -g -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -fno-omit-frame-pointer")


# ──────────────────────────────────────────────────────────────────
# 1. Optional sanitiser toggle
# ──────────────────────────────────────────────────────────────────
option(ENABLE_SANITIZERS "Build with Address + UB Sanitiser" OFF)
if(ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        -fsanitize=address,undefined
        -fno-pie                   # ASan shadow-range fix
    )
    add_link_options(
        -fsanitize=address,undefined
        -no-pie
    )
endif()

# ──────────────────────────────────────────────────────────────────
# 2. 3rd-party libraries
# ──────────────────────────────────────────────────────────────────
include(FetchContent)

# 2.1 spdlog (header-only)
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG        v1.15.3
)
FetchContent_MakeAvailable(spdlog)

# 2.2 OGDF (pre-built – point to your prefix)
set(OGDF_DIR "/home/politov/OGDF" CACHE PATH "OGDF install prefix")
find_package(OGDF CONFIG REQUIRED PATHS "${OGDF_DIR}/lib/cmake")

# 2.3 OpenMP
find_package(OpenMP REQUIRED)

# ──────────────────────────────────────────────────────────────────
# 3. Project libraries
# ──────────────────────────────────────────────────────────────────
add_subdirectory(src/util)   # → libsbutil
add_subdirectory(src/io)     # → libsbio

# ──────────────────────────────────────────────────────────────────
# 4. Main executable
# ──────────────────────────────────────────────────────────────────
add_executable(sbfind src/main.cpp)
target_include_directories(sbfind PRIVATE src)
target_compile_options(sbfind PRIVATE -O2  -g -fno-omit-frame-pointer)


target_link_libraries(sbfind PRIVATE
    sbio
    sbutil
    OGDF
    spdlog::spdlog_header_only
    OpenMP::OpenMP_CXX
)

# ──────────────────────────────────────────────────────────────────
# 5. Warnings & optimisation for *this* target only
# ──────────────────────────────────────────────────────────────────
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(sbfind PRIVATE
        -Wall -Wextra -Wshadow -pedantic
    )
endif()

                                
                                
                                # cmake_minimum_required(VERSION 3.20)
                                # project(sbSPQRnew LANGUAGES CXX)

                                # # ------------------------------------------------------------------#
                                # # Global options                                                    #
                                # # ------------------------------------------------------------------#
                                # set(CMAKE_CXX_STANDARD 20)
                                # set(CMAKE_CXX_STANDARD_REQUIRED ON)

                                # option(ENABLE_SANITIZERS
                                #     "Build with AddressSanitizer + UndefinedBehaviourSanitizer" OFF)




                                # # ------------------------------------------------------------------#
                                # # Sanitizer flags (GCC / Clang)                                     #
                                # # ------------------------------------------------------------------#
                                # if(ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
                                #     add_compile_options(
                                #         -fsanitize=address
                                #         -fsanitize=undefined
                                #         -fno-omit-frame-pointer    # keep full back-traces
                                #         -fno-pie                   # avoid ASan shadow-range clash
                                #         -g                         # debug symbols
                                #     )
                                #     add_link_options(
                                #         -fsanitize=address
                                #         -fsanitize=undefined
                                #         -no-pie                    # link a non-PIE executable
                                #     )
                                #     message(STATUS "Building with AddressSanitizer + UBSan (no-PIE)")
                                # endif()

                                # # ------------------------------------------------------------------#
                                # # 3rd-party: spdlog (header-only)                                   #
                                # # ------------------------------------------------------------------#
                                # include(FetchContent)
                                # FetchContent_Declare(
                                # spdlog
                                # GIT_REPOSITORY https://github.com/gabime/spdlog.git
                                # GIT_TAG        v1.15.3
                                # )
                                # FetchContent_MakeAvailable(spdlog)

                                # # ------------------------------------------------------------------#
                                # # 3rd-party: OGDF (built/installed elsewhere)                       #
                                # # ------------------------------------------------------------------#
                                # set(OGDF_DIR "/home/politov/OGDF" CACHE PATH "OGDF install prefix")
                                # find_package(OGDF CONFIG REQUIRED
                                #             PATHS "${OGDF_DIR}/lib/cmake")

                                # find_package(OpenMP REQUIRED)


                                # # ------------------------------------------------------------------#
                                # # Project subdirectories (static libs)                              #
                                # # ------------------------------------------------------------------#
                                # add_subdirectory(src/util)   # → libsbutil
                                # add_subdirectory(src/io)     # → libsbio

                                # # ------------------------------------------------------------------#
                                # # Main executable                                                   #
                                # # ------------------------------------------------------------------#
                                # add_executable(sbfind src/main.cpp)
                                # target_include_directories(sbfind PRIVATE src)   # #include "util/…"
                                # target_link_libraries(sbfind PRIVATE
                                #     sbio
                                #     sbutil
                                #     OGDF
                                #     spdlog::spdlog_header_only
                                #     OpenMP::OpenMP_CXX
                                # )

                                # # target_compile_options(sbfind PRIVATE -fopenmp)
                                # # target_link_libraries(sbfind PRIVATE -fopenmp)


                                # # ------------------------------------------------------------------#
                                # # Nice warnings                                                     #
                                # # ------------------------------------------------------------------#
                                # if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
                                #     # target_compile_options(sbfind PRIVATE -O3 -Wall -Wextra -Wshadow -pedantic -fopenmp -g -march=native -flto -fno-omit-frame-pointer)
                                #     target_compile_options(sbfind PRIVATE -O2 -g -fno-omit-frame-pointer -DNDEBUG)

                                # endif()
