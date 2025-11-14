find_package(SDL2 QUIET)
find_package(SDL2_image QUIET)
find_package(SDL2_ttf QUIET)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(COMPILE_SDL2_CODE TRUE)
if(${SDL2_ttf_FOUND})
    message("SDL2_ttf found")
else()
    message("SDL2_ttf NOT found!")
    set(COMPILE_SDL2_CODE FALSE)
endif()

if(${SDL2_image_FOUND})
    message("SDL2_image found")
else()
    message("SDL2_image NOT found!")
    set(COMPILE_SDL2_CODE FALSE)
endif()

set(COMMON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/json_config.cpp
                   ${CMAKE_CURRENT_SOURCE_DIR}/src/common.cpp)

if(${COMPILE_SDL2_CODE})
    add_compile_definitions(WITH_SDL2)

    include_directories(${SDL2_INCLUDE_DIRS})

    set(COMMON_SOURCES ${COMMON_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/SDL_plot.cpp)
    set(COMMON_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/include/SDL_plot.hpp)

    set(TOOLS_LIBRARIES ${SDL2_LIBRARIES} SDL2_image::SDL2_image SDL2_ttf::SDL2_ttf)
endif()

find_package(nlohmann_json 3.2.0 REQUIRED)
find_package( Boost COMPONENTS program_options REQUIRED )

set(TOOLS_LIBRARIES ${TOOLS_LIBRARIES} ${BUILD_LIB_TARGET} ${Boost_LIBRARIES})

set(SPECIES_SIM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/species_sim.cpp
                        ${COMMON_SOURCES})

set(TISSUE_SAMPLER ${CMAKE_CURRENT_SOURCE_DIR}/src/tissue_sampler.cpp
                    ${COMMON_SOURCES})

set(DESCENDANTS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/descendants_builder.cpp
                ${COMMON_SOURCES})

set(MUTATIONS_SIM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/mutations_sim.cpp
                        ${COMMON_SOURCES})

set(BUILD_SBS_CONTEXT_INDEX_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/build_SBS_context_index.cpp
                                ${COMMON_SOURCES})

set(BUILD_ID_CONTEXT_INDEX_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/build_ID_context_index.cpp
                        ${COMMON_SOURCES})

set(SAMPLE_INDEX_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/sample_context_index.cpp)

add_executable(species_sim ${SPECIES_SIM_SOURCES})
target_include_directories(species_sim PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(species_sim ${TOOLS_LIBRARIES} nlohmann_json::nlohmann_json)
install(TARGETS species_sim)

add_executable(descendants_builder ${DESCENDANTS_SOURCES})
target_include_directories(descendants_builder PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(descendants_builder ${TOOLS_LIBRARIES})
install(TARGETS descendants_builder)

add_executable(tissue_sampler ${TISSUE_SAMPLER})
target_include_directories(tissue_sampler PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(tissue_sampler ${TOOLS_LIBRARIES} nlohmann_json::nlohmann_json)
install(TARGETS tissue_sampler)

add_executable(build_SBS_context_index ${BUILD_SBS_CONTEXT_INDEX_SOURCES})
target_include_directories(build_SBS_context_index PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(build_SBS_context_index ${TOOLS_LIBRARIES})
install(TARGETS build_SBS_context_index)

add_executable(build_ID_context_index ${BUILD_ID_CONTEXT_INDEX_SOURCES})
target_include_directories(build_ID_context_index PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(build_ID_context_index ${TOOLS_LIBRARIES})
install(TARGETS build_ID_context_index)

add_executable(mutations_sim ${MUTATIONS_SIM_SOURCES})
target_include_directories(mutations_sim PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(mutations_sim ${TOOLS_LIBRARIES} nlohmann_json::nlohmann_json)
install(TARGETS mutations_sim)