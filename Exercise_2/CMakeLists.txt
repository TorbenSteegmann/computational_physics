cmake_minimum_required(VERSION 3.5)

project(exercise_2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(exercise_2
   main.cc
   ../Matplotlib/matplotlibcpp.h
)

target_link_libraries(exercise_2 PRIVATE "path\\to\\pythonxxx.lib")

target_include_directories(exercise_2 PRIVATE "path\\to\\Pythonxxx\\include")



include(GNUInstallDirs)
install(TARGETS exercise_2
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
