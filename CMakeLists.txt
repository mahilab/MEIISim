cmake_minimum_required(VERSION 3.0.0)
project(OpenWristSim VERSION 0.1.0)

include(FetchContent) 
FetchContent_Declare(mahi-util GIT_REPOSITORY https://github.com/mahilab/mahi-util.git) 
FetchContent_MakeAvailable(mahi-util)
FetchContent_Declare(mahi-robo GIT_REPOSITORY https://github.com/mahilab/mahi-robo.git) 
FetchContent_MakeAvailable(mahi-robo)

include_directories(include)

add_executable(test src/test.cpp)
target_link_libraries(test mahi::util mahi::robo)
