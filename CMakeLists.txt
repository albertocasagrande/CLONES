cmake_minimum_required(VERSION 3.23)

project(RACES VERSION 1.1.0
              DESCRIPTION "RACES is an Advanced Cancer Evolution Simulation"
              HOMEPAGE_URL https://github.com/albertocasagrande/RACES
              LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_MACOSX_RPATH 1)
endif()
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.30.0")
cmake_policy(SET CMP0167 NEW)
endif()
add_compile_options(-Wall -Wextra -Wpedantic)

set(RACES_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/allele.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/archive.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/basic_IO.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/binary_logger.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/cell.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/cell_event.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/cna.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/csv_reader.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/descendant_forest.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/driver_storage.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/ending_conditions.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/event_wrapper.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/fasta_chr_reader.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/fasta_utils.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/genome_mutations.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/genomic_position.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/genomic_region.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/genomic_sequence.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/germline.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/id_context.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/id_context_index.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/id_signature.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/lineage_graph.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/logger.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/logics.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/mutant_mutation.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/mutant_properties.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/mutation_engine.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/mutation_list.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/mutation.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/mutational_properties.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/palette.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/phyloXML.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/phylogenetic_forest.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/plot_2D.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/position.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/position_set.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/progress_bar.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/rate_update.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/read.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/read_simulator.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/sample_specification.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/sampling.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/sbs_context.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/sbs_signature.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/sequencer.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/sid.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation_event.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/species.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/statistics.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/timed_event.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/tissue.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/tissue_sample.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/utils.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/warning.cpp
)

set(RACES_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/include/allele.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/archive.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/basic_IO.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/binary_logger.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/bucket.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/cell.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/cell_event.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/cna.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/csv_reader.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/descendant_forest.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/driver_storage.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/ending_conditions.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/event_wrapper.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/filter.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/fasta_reader.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/fasta_chr_reader.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/fasta_utils.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/genome_mutations.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/genomic_position.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/genomic_region.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/genomic_sequence.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/germline.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/id_context.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/id_context_index.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/id_signature.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/imap.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/index.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/irb_tree.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/label_tour.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/lineage_graph.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/logger.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/logics.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/mutant_id.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/mutant_mutation.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/mutant_properties.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/mutation_engine.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/mutation_list.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/mutation_spec.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/mutation.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/mutational_properties.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/palette.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/palette.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/phyloXML.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/phylogenetic_forest.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/plot_2D.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/position.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/position_set.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/progress_bar.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/rate_update.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/read.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/read_simulator.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/sample_specification.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/sampling.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/sbs_context.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/sbs_context_index.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/sbs_signature.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/sequencer.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/sid.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/signature.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/simulation.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/simulation_event.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/species.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/statistics.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/time.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/timed_event.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/tissue.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/tissue_plotter.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/tissue_sample.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/variables.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/utils.hpp
                ${CMAKE_CURRENT_SOURCE_DIR}/include/warning.hpp
)

set(WITH_MATPLOT false)

if(NOT DEFINED USE_OPTIONAL_LIBS)
    set(USE_OPTIONAL_LIBS TRUE)
endif()

if(${USE_OPTIONAL_LIBS})
    find_package(Matplot++ QUIET)

    if(${Matplot++_FOUND})
        message("Matplot++ found")

        set(WITH_MATPLOT true)

        set(RACES_LIBRARIES ${RACES_LIBRARIES} Matplot++::matplot)
    else()
        message("Matplot++ NOT found")
    endif()
else()
    message("Avoiding optional libraries")
endif()

set(WITH_INDICATORS false)

if(NOT EXISTS ${CMAKE_SOURCE_DIR}/include/indicators.hpp)
    message("Downloading \"p-ranav/indicators\" library...")
    file(DOWNLOAD
        https://raw.githubusercontent.com/p-ranav/indicators/master/single_include/indicators/indicators.hpp
        ${CMAKE_SOURCE_DIR}/include/indicators.hpp
        STATUS DOWNLOAD_RESULT)
    list(GET DOWNLOAD_RESULT 1 DOWNLOAD_RESULT_MESSAGE)
    if(${DOWNLOAD_RESULT_MESSAGE} STREQUAL "\"No error\"")
        set(WITH_INDICATORS true)
        message("done")
        message("Using \"p-ranav/indicators\"")

        set(RACES_HEADER ${RACES_HEADER} ${CMAKE_CURRENT_SOURCE_DIR}/include/indicators.hpp)
    else()
        set(WITH_INDICATORS false)
        message("Download error: ${DOWNLOAD_RESULT_MESSAGE}")
        message("Avoiding \"p-ranav/indicators\"")
    endif()
else()
    set(WITH_INDICATORS true)
    message("Using \"p-ranav/indicators\"")

    set(RACES_HEADER ${RACES_HEADER} ${CMAKE_CURRENT_SOURCE_DIR}/include/indicators.hpp)
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/variables.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/include/variables.hpp @ONLY)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_library(objlib OBJECT ${RACES_SOURCES})
target_link_libraries(objlib ${RACES_LIBRARIES})

set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

if(NOT DEFINED BUILD_STATIC_LIBS)
    set(BUILD_STATIC_LIBS TRUE)
endif()

if(${BUILD_STATIC_LIBS})
    set(BUILD_LIB_TARGET libRACES_static)
    set(LIB_INSTALL_TARGETS libRACES_static)
endif()

if(NOT DEFINED BUILD_ONLY_LIBRARIES)
    set(BUILD_ONLY_LIBRARIES FALSE)
endif()

if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_LIB_TARGET libRACES)
    set(BUILD_SHARED_LIBS TRUE)
endif()

if(${BUILD_SHARED_LIBS})
    set(LIB_INSTALL_TARGETS libRACES ${LIB_INSTALL_TARGETS})
endif()

if((NOT ${BUILD_SHARED_LIBS}) AND (NOT ${BUILD_STATIC_LIBS}))
    message( FATAL_ERROR "Either static or dynamic library must be built." )
endif()

if(NOT ${BUILD_ONLY_LIBRARIES})
    message("Building executable, Python wrapper, and tests")

    add_subdirectory(python)
    add_subdirectory(tools)

    enable_testing()

    add_subdirectory(tests)
else()
    message("Building exclusively RACES libraries")
endif()

if(${BUILD_SHARED_LIBS})
    add_library(libRACES SHARED $<TARGET_OBJECTS:objlib>)
    target_link_libraries(libRACES ${RACES_LIBRARIES})
    set_target_properties(libRACES PROPERTIES VERSION ${PROJECT_VERSION} OUTPUT_NAME RACES)
endif()

if(${BUILD_STATIC_LIBS})
    add_library(libRACES_static STATIC $<TARGET_OBJECTS:objlib>)
    target_link_libraries(libRACES_static ${RACES_LIBRARIES})
    set_target_properties(libRACES_static PROPERTIES VERSION ${PROJECT_VERSION} OUTPUT_NAME RACES)
endif()

if(NOT DEFINED INSTALL_HEADERS)
    set(INSTALL_HEADERS TRUE)
endif()

if(${INSTALL_HEADERS})
  if(${BUILD_SHARED_LIBS})
    set_property(TARGET libRACES PROPERTY PUBLIC_HEADER ${RACES_HEADER})
  endif()

  if(${BUILD_STATIC_LIBS})
     set_property(TARGET libRACES_static PROPERTY PUBLIC_HEADER ${RACES_HEADER})
  endif()
  install(TARGETS ${LIB_INSTALL_TARGETS} FILE_SET HEADERS)
else()
  install(TARGETS ${LIB_INSTALL_TARGETS})
endif()
