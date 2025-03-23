cmake_minimum_required(VERSION 3.10)
project(MTR_Detect VERSION 0.1 DESCRIPTION "A Hello World Project")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file (GLOB_RECURSE SOURCE_FILES src/*.cpp)

# 将 pybind11 添加为子目录
add_subdirectory(lib/pybind11)

# 查找 Python
find_package(Python3 COMPONENTS Development NumPy REQUIRED)

# This project will output an executable file
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# 添加pthread链接
target_link_libraries(${PROJECT_NAME} PRIVATE pthread)

# Create a simple configuration header
configure_file(config.h.in config.h)

# Include the configuration header in the build
target_include_directories(${PROJECT_NAME}
    PRIVATE
        "${PROJECT_BINARY_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        ${Python3_INCLUDE_DIRS}
        ${Python3_NumPy_INCLUDE_DIRS}
        "${CMAKE_SOURCE_DIR}/third_party/pybind11/include"
        "/home/rose/soft/matplotlib-cpp/"
)


target_link_libraries(${PROJECT_NAME}
    PRIVATE
        pybind11::module
        ${Python3_LIBRARIES}
        ${Python3_NumPy_LIBRARIES}
)