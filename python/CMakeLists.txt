find_package(Python COMPONENTS Interpreter Development)
find_package(Boost COMPONENTS python)

set(COMPILE_PYTHON_BINDINGS TRUE)

if(NOT ${Boost_python_FOUND})
message("Boost.Python is missing: Python bindings will NOT be compiled")
set(COMPILE_PYTHON_BINDINGS FALSE)
endif()

if(NOT ${PYTHONLIBS_FOUND})
message("Python libraries are missing: Python bindings will NOT be compiled")
set(COMPILE_PYTHON_BINDINGS FALSE)
endif()

if(${COMPILE_PYTHON_BINDINGS})
    message("Python bindings will be compiled")

    # include directories
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
    include_directories(Python::Module)
    include_directories(${Boost_INCLUDE_DIRS})

    # create the lib
    add_library(pyCLONES MODULE
                      ${CMAKE_CURRENT_SOURCE_DIR}/src/bindings.cpp
                      ${CMAKE_CURRENT_SOURCE_DIR}/src/mutant.cpp  
                      ${CMAKE_CURRENT_SOURCE_DIR}/src/epigenetic_rates.cpp  
                      ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation_wrapper.cpp)

    set_target_properties(pyCLONES PROPERTIES OUTPUT_NAME CLONES PREFIX "")

    # link
    target_link_libraries(pyCLONES ${Boost_LIBRARIES} Python::Module ${BUILD_LIB_TARGET})

    execute_process(COMMAND ${Python_EXECUTABLE} -m pytest --version
        OUTPUT_VARIABLE PYTEST_output
        ERROR_VARIABLE  PYTEST_error
        RESULT_VARIABLE PYTEST_result)
    if(NOT ${PYTEST_result} EQUAL 0)
        message("python-test will not be executed: pytest package not available")
    else()
        enable_testing()
        add_test (NAME python-tests
            COMMAND Python::Interpreter -m pytest # Or just COMMAND pytest
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/
            )
        set_tests_properties(python-tests
            PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}:$ENV{PYTHONPATH}")
    endif()
endif()
