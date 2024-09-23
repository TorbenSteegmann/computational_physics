cmake_minimum_required(VERSION 3.5)

project(exam LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(exam
   main.cc
   ../Matplotlib/matplotlibcpp.h
)

target_link_libraries(exam PRIVATE "path\\to\\pythonxxx.lib")

target_include_directories(exam PRIVATE "path\\to\\Pythonxxx\\include")

include(GNUInstallDirs)
install(TARGETS exam
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
