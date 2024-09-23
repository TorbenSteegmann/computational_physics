cmake_minimum_required(VERSION 3.5)

project(exercise_5 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(exercise_5
   main.cc
   ../Matplotlib/matplotlibcpp.h
)

target_link_libraries(exercise_5 PRIVATE "path\\to\\pythonxxx.lib")

target_include_directories(exercise_5 PRIVATE "path\\to\\Pythonxxx\\include")

include(GNUInstallDirs)
install(TARGETS exercise_5
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
