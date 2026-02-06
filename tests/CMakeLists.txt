find_package(Boost COMPONENTS unit_test_framework)

if(${Boost_unit_test_framework_FOUND})

    message("C++ code will tested")

    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    set(TEST_LIBRARIES ${Boost_LIBRARIES} ${BUILD_LIB_TARGET} )

    set(CLONES_TESTS archive fasta_reader context_index genome_mutations)

    set(ASSET_DIR "${CMAKE_CURRENT_SOURCE_DIR}/assets")
    set(SBS_EXAMPLE "${ASSET_DIR}/v3.3_SBS3_PROFILE.txt")
    if(NOT EXISTS ${SBS_EXAMPLE})
        message("Downloading SBS profile example from sanger...")
        file(DOWNLOAD
            https://cancer.sanger.ac.uk/signatures/documents/1944/v3.3_SBS3_PROFILE.txt
            ${SBS_EXAMPLE}
            STATUS DOWNLOAD_RESULT)

        list(GET DOWNLOAD_RESULT 1 DOWNLOAD_RESULT_MESSAGE)
        if(${DOWNLOAD_RESULT_MESSAGE} STREQUAL "\"No error\"")
            add_executable( test_sbs_signature sbs_signature.cpp)
            add_dependencies( test_sbs_signature ${BUILD_LIB_TARGET} )
            target_link_libraries(test_sbs_signature ${TEST_LIBRARIES})
            target_compile_definitions(test_sbs_signature PUBLIC SBS_EXAMPLE="${SBS_EXAMPLE}")
        endif()
    else()
        add_executable( test_sbs_signature sbs_signature.cpp)
        add_dependencies( test_sbs_signature ${BUILD_LIB_TARGET} )
        target_link_libraries(test_sbs_signature ${TEST_LIBRARIES})
        target_compile_definitions(test_sbs_signature PUBLIC SBS_EXAMPLE="${SBS_EXAMPLE}")
    endif()
    set(ID_EXAMPLE "${ASSET_DIR}/COSMIC_v3.4_ID_GRCh37.txt")
    if(NOT EXISTS ${ID_EXAMPLE})
        message("Downloading ID profile example from sanger...")
        file(DOWNLOAD
        https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.4_ID_GRCh37.txt
            ${ID_EXAMPLE}
            STATUS DOWNLOAD_RESULT)

        list(GET DOWNLOAD_RESULT 1 DOWNLOAD_RESULT_MESSAGE)
        if(${DOWNLOAD_RESULT_MESSAGE} STREQUAL "\"No error\"")
            add_executable( test_id_signature id_signature.cpp)
            add_dependencies( test_id_signature ${BUILD_LIB_TARGET} )
            target_link_libraries(test_id_signature ${TEST_LIBRARIES})
            target_compile_definitions(test_id_signature PUBLIC ID_EXAMPLE="${ID_EXAMPLE}")
        endif()
    else()
        add_executable( test_id_signature id_signature.cpp )
        add_dependencies( test_id_signature ${BUILD_LIB_TARGET} )
        target_link_libraries(test_id_signature ${TEST_LIBRARIES})
        target_compile_definitions(test_id_signature PUBLIC ID_EXAMPLE="${ID_EXAMPLE}")
    endif()

    foreach(TEST ${CLONES_TESTS})
        add_executable( test_${TEST} ${TEST}.cpp)
        add_dependencies( test_${TEST} ${BUILD_LIB_TARGET} )
        target_link_libraries(test_${TEST} ${TEST_LIBRARIES})
        ADD_TEST(${TEST} test_${TEST})
    endforeach()

    target_compile_definitions(test_fasta_reader PUBLIC FASTA_FILE="${ASSET_DIR}/fasta_test.fasta"
    							                        FASTA_INDEX_ERR="${ASSET_DIR}/fasta_test_index_err.fasta")
    target_compile_definitions(test_context_index PUBLIC FASTA_FILE="${ASSET_DIR}/context_index_test.fasta")
    target_compile_definitions(test_genome_mutations PUBLIC FASTA_FILE="${ASSET_DIR}/context_index_test.fasta")
endif()
