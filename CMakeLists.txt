cmake_minimum_required(VERSION 3.2.0)
# With CMAKE >= 3.2 external projects support byproducts, this is needed
# to support Ninja as generator

project(BioDynaMo)
set(BDM_VERSION 0.1.0)

# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(test "Build all tests." ON) # Makes boolean 'test' available.
if(APPLE)
   # valgrind currently broken on macOS 10.12
   option(valgrind "Enable valgrind tests and make build compatible with valgrind tool." OFF)
else()
   option(valgrind "Enable valgrind tests and make build compatible with valgrind tool." ON)
endif()
option(coverage "Enable test coverage report generation. Sets build type to coverage" OFF)

set(BDM_CMAKE_DIR cmake)

if(UNIX AND NOT APPLE)
  set(LINUX TRUE)
endif()

add_definitions(-DBDM_SRC_DIR=\"${CMAKE_SOURCE_DIR}/src\")  # escaped quotes are mandatory

# -------------------- find packages ------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
# CMake cannot find LLVM's Clang OpenMP library by default (this makes
# find_package(OpenMP) fail). Therefore we manually specify the OpenMP flags.
# Furthermore, rootcling cannot find omp.h by default, so we copy this (with a
# shell script to /usr/local/Cellar/biodynamo, which is brew's default install
# directory), and include this directory. We cannot directly include the
# original directory, because of header conflicts (such as stdint.h)
if(APPLE)
  link_directories(/usr/local/opt/llvm/lib)
  set(OpenMP_C_FLAGS -fopenmp=libomp)
  set(OpenMP_CXX_FLAGS -fopenmp=libomp)
  include_directories(/usr/local/Cellar/biodynamo)
else()
  find_package(OpenMP REQUIRED)
endif()
if (test)
  # must be before include(${ROOT_USE_FILE}) otherwise test compile fails
  find_package(Threads REQUIRED)
endif()

find_package(Git)
find_package(ROOT REQUIRED COMPONENTS Eve Geom Gui)

find_package(ClangTools)
if ("$ENV{CMAKE_EXPORT_COMPILE_COMMANDS}" STREQUAL "1" OR CLANG_TIDY_FOUND)
  # Generate a Clang compile_commands.json "compilation database" file for use
  # with various development tools, such as Vim's YouCompleteMe plugin.
  # See http://clang.llvm.org/docs/JSONCompilationDatabase.html
  set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
endif()

find_package(ParaView 5.3 QUIET OPTIONAL_COMPONENTS vtkPVPythonCatalyst
  vtkPVVTKExtensionsDefault vtkPVServerManagerRendering)
if(${ParaView_FOUND})
  message(STATUS "Found ParaView")
  include("${PARAVIEW_USE_FILE}")
  add_definitions("-DUSE_CATALYST")
else()
  message(WARNING "ParaView not found")
endif()

find_package(VTune)
if(${VTune_FOUND})
  include_directories(${VTUNE_INCLUDE_DIR})
  add_definitions("-DUSE_VTUNE")
else()
  message(WARNING "VTune not found")
endif()

# -------------------- includes -----------------------------------------------
include(BioDynaMo)
include(ExternalProject)
include(cmake/doxygen.cmake)

# -------------------- set default build type and compiler flags ---------------
if(coverage)
  set(CMAKE_BUILD_TYPE Coverage)
elseif(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
include(cmake/SetCompilerFlags.cmake)

# -------------------- set rpath options ---------------------------------------
# When building, use the RPATH
set(CMAKE_SKIP_BUILD_RPATH FALSE)           # don't skip the full RPATH for the build tree
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
if(APPLE)
  set(CMAKE_INSTALL_RPATH
	"/usr/local/lib/paraview-5.4"
	"/usr/local/opt/qt5/lib"
	"/usr/local/lib/biodynamo")
endif()

# -------------------- setup mpark variant library -----------------------------
set(MPARK_LIB "${CMAKE_CURRENT_BINARY_DIR}/mpark/")
file(MAKE_DIRECTORY "${MPARK_LIB}/")
execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf "${CMAKE_SOURCE_DIR}/third_party/mpark_4ee0086.zip"
                WORKING_DIRECTORY "${MPARK_LIB}")

# -------------------- build shared library libbiodynamo -----------------------
# add the headers into the project
include(${ROOT_USE_FILE})
include_directories("${CMAKE_SOURCE_DIR}/src")
include_directories("${CMAKE_SOURCE_DIR}/test")
include_directories("${CMAKE_SOURCE_DIR}/third_party")
include_directories("${CMAKE_SOURCE_DIR}/demo")
include_directories("${MPARK_LIB}")

# build libbiodynamo
file(GLOB_RECURSE HEADERS "${CMAKE_SOURCE_DIR}/src/*.h")
file(GLOB_RECURSE LIB_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cc")
build_libbiodynamo(biodynamo
                   SOURCES ${LIB_SOURCES}
                   HEADERS ${HEADERS}
                   LIBRARIES ${ROOT_LIBRARIES})
if(${VTune_FOUND})
  target_link_libraries(biodynamo ${VTUNE_LIBRARIES})
endif()

# add custom clean target for shared library
add_custom_target(libbdmclean COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/CMakeFiles/biodynamo.dir/cmake_clean.cmake")
# custom clean target for biodynamo that excludes ExternalProjects
add_custom_target(bdmclean)
add_dependencies(bdmclean libbdmclean)

generate_rootlogon()
if (${ParaView_FOUND})
  link_directories($ENV{ParaView_DIR}/../../paraview-5.4)
  target_link_libraries(biodynamo vtkPVPythonCatalyst vtkIOParallelXML
    vtkPVVTKExtensionsDefault vtkPVServerManagerRendering)
endif()

# create cmake files used by external builds
configure_file("cmake/UseBioDynaMo.cmake.in" "UseBioDynaMo.cmake" @ONLY)


# -------------------- build test project  ------------------------------------
if (test)
  enable_testing()

  include(cmake/Testing.cmake)

  # unit tests
  #   main test executable
  file(GLOB TEST_SOURCES ${CMAKE_SOURCE_DIR}/test/unit/*.cc
                         ${CMAKE_SOURCE_DIR}/test/unit/biology_module/*.cc)
  file(GLOB TEST_HEADERS ${CMAKE_SOURCE_DIR}/test/unit/*.h
                         ${CMAKE_SOURCE_DIR}/test/unit/biology_module/*.h)
  bdm_add_test_executable(runBiodynamoTestsMain
                          SOURCES ${TEST_SOURCES}
                          HEADERS ${TEST_HEADERS})

  #   separate test binaries - each source file is compiled into a separate binary
  file(GLOB TEST_SOURCES ${CMAKE_SOURCE_DIR}/test/unit/separate_binary/*.cc)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/empty_header.h "")
  foreach(file ${TEST_SOURCES})
    get_filename_component(filename ${file} NAME_WE)
    set(TEST_HEADER ${CMAKE_SOURCE_DIR}/test/unit/separate_binary/${filename}.h)
    if (NOT EXISTS ${TEST_HEADER})
      set(TEST_HEADER ${CMAKE_CURRENT_BINARY_DIR}/empty_header.h)
    endif()
    bdm_add_test_executable("runBiodynamoTests_${filename}"
                            SOURCES ${file}
                            HEADERS ${TEST_HEADER})
  endforeach()

  # integration tests
  if (NOT coverage)
    #   backup-restore
    bdm_add_executable(backup-restore
                      SOURCES test/integration/backup_restore.cc
                      HEADERS test/integration/backup_restore.h
                      LIBRARIES biodynamo)
    add_dependencies(check backup-restore)
    add_test(NAME "integration-backup-restore"
             COMMAND ${CMAKE_SOURCE_DIR}/test/integration/backup_restore.sh)
    bdm_add_executable(soma_clustering
                      SOURCES test/integration/soma_clustering.cc
                      HEADERS test/integration/soma_clustering.h
                      LIBRARIES biodynamo)
    add_dependencies(check soma_clustering)
    add_test(NAME "integration-soma-clustering"
             COMMAND ${CMAKE_SOURCE_DIR}/build/soma_clustering)
  endif()
endif()

# -------------------- Install configuration -----------------------------------
include(cmake/Installation.cmake)

# -------------------- add targets for code style and submission checks --------
include(cmake/CppStyleGuideChecks.cmake)

# -------------------- build snap package --------------------------------------
# create snapcraft config file
configure_file("cmake/snapcraft.yaml.in" "snapcraft.yaml" @ONLY)

add_custom_target(snap-package
  COMMAND ../cmake/build_snap_package.sh ${CMAKE_INSTALL_PREFIX}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Build snap package")

# -------------------- build cell_division_module executable -------------------
bdm_add_executable(cell_division_module
                   SOURCES demo/cell_division_module.cc
                   HEADERS demo/cell_division_module.h
                   LIBRARIES biodynamo)

# -------------------- build cell_diffusion_module executable ------------------
bdm_add_executable(diffusion_module
                   SOURCES demo/diffusion_module.cc
                   HEADERS demo/diffusion_module.h
                   LIBRARIES biodynamo)

# --------------------- build cell_diffusion_module executable ------------------
bdm_add_executable(cell_gene_module
                   SOURCES demo/cell_gene_module.cc
                   HEADERS demo/cell_gene_module.h
                   LIBRARIES biodynamo)
