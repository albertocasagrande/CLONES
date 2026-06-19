find_package(Boost COMPONENTS unit_test_framework)

function(download_asset_and_add_test test_name asset_filename asset_url asset_variable)
    set(FULL_TEST_NAME "test_${test_name}")
    if(NOT EXISTS ${asset_filename})
        message("Downloading assets for ${test_name}...")
        file(DOWNLOAD
            ${asset_url}
            ${asset_filename}
            STATUS DOWNLOAD_RESULT)

        list(GET DOWNLOAD_RESULT 1 DOWNLOAD_RESULT_MESSAGE)
        if(${DOWNLOAD_RESULT_MESSAGE} STREQUAL "\"No error\"")

            add_executable( ${FULL_TEST_NAME} ${test_name}.cpp)
            add_dependencies( ${FULL_TEST_NAME} ${BUILD_LIB_TARGET} )
            target_link_libraries(${FULL_TEST_NAME} ${TEST_LIBRARIES})
            target_compile_definitions(${FULL_TEST_NAME} PUBLIC ${asset_variable}="${asset_filename}")
            ADD_TEST(${test_name} ${FULL_TEST_NAME})
        endif()
    else()
        add_executable( ${FULL_TEST_NAME} ${test_name}.cpp)
        add_dependencies( ${FULL_TEST_NAME} ${BUILD_LIB_TARGET} )
        target_link_libraries(${FULL_TEST_NAME} ${TEST_LIBRARIES})
        target_compile_definitions(${FULL_TEST_NAME} PUBLIC ${asset_variable}="${asset_filename}")
        ADD_TEST(${test_name} ${FULL_TEST_NAME})
    endif()
endfunction()

if(${Boost_unit_test_framework_FOUND})

    message("C++ code will tested")

    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    set(TEST_LIBRARIES ${Boost_LIBRARIES} ${BUILD_LIB_TARGET} )

    set(RACES_TESTS archive fasta_reader context_index genome_mutations bucket index tissue_simulation)

    foreach(TEST ${RACES_TESTS})
        add_executable( test_${TEST} ${TEST}.cpp test_common.cpp)
        add_dependencies( test_${TEST} ${BUILD_LIB_TARGET} )
        target_link_libraries(test_${TEST} ${TEST_LIBRARIES})
        ADD_TEST(${TEST} test_${TEST})
    endforeach()

    set(ASSET_DIR "${CMAKE_CURRENT_SOURCE_DIR}/assets")
    target_compile_definitions(test_fasta_reader PUBLIC FASTA_FILE="${ASSET_DIR}/fasta_test.fasta"
    							                        FASTA_INDEX_ERR="${ASSET_DIR}/fasta_test_index_err.fasta")
    target_compile_definitions(test_context_index PUBLIC FASTA_FILE="${ASSET_DIR}/context_index_test.fasta")
    target_compile_definitions(test_genome_mutations PUBLIC FASTA_FILE="${ASSET_DIR}/context_index_test.fasta")

    download_asset_and_add_test(sbs_signature "${ASSET_DIR}/v3.3_SBS3_PROFILE.txt" 
                                "https://cancer.sanger.ac.uk/signatures/documents/1944/v3.3_SBS3_PROFILE.txt"
                                "SBS_EXAMPLE")
    download_asset_and_add_test(id_signature "${ASSET_DIR}/COSMIC_v3.4_ID_GRCh37.txt" 
                                "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.4_ID_GRCh37.txt"
                                "ID_EXAMPLE")

endif()
