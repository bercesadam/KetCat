cmake_minimum_required(VERSION 3.16)
project(MultiExeProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ----------------------------
# Collect all cpp files
# ----------------------------
file(GLOB_RECURSE ALL_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)

# ----------------------------
# Compiler-specific settings
# ----------------------------
if (MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -g)
endif()

if (MSVC)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS _CRT_NONSTDC_NO_DEPRECATE)
endif()

# ----------------------------
# Create one executable per .cpp
# ----------------------------
foreach(SRC_FILE ${ALL_SOURCES})

    file(RELATIVE_PATH REL_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${SRC_FILE}
    )

    string(REPLACE "/" "_" EXE_NAME ${REL_PATH})
    string(REPLACE ".cpp" "" EXE_NAME ${EXE_NAME})

    add_executable(${EXE_NAME} ${SRC_FILE} "include/compiler/physical_instruction.h" "include/quantum_processor/simulation_observer.h")

    if (MSVC)
        set_target_properties(${EXE_NAME} PROPERTIES LINK_FLAGS "/STACK:67108864")
    else()
        # Linux/GCC
        target_link_options(${EXE_NAME} PRIVATE "-Wl,-z,stack-size=67108864")
    endif()

    target_include_directories(${EXE_NAME}
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )

    target_link_libraries(${EXE_NAME}
        PRIVATE
    )

endforeach()
