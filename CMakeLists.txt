project(Analysis)
cmake_minimum_required(VERSION 2.8)

#ROOT libraries
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO  MINUIT)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

#Generate the dictionary for vector<vector>
ROOT_GENERATE_DICTIONARY( G__Vector LINKDEF ../inc/LinkDef.h )

add_library(Vector SHARED G__Vector.cxx)
target_link_libraries(Vector ${ROOT_LIBRARIES})

#Bring the headers
include_directories(inc)
include_directories($ENV{EIGEN_INSTALL})

#Add the CMakeLists.txt on the subdirectory
add_subdirectory(src ${PROJECT_SOURCE_DIR}/build/sourcebuild)

#Compile custom libraries
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/inc/*.h)

# COMPIRE CUSTOM SOURCES
#IN FOLDER WAVEFORM ############################################################
add_executable(waveformAverage ${PROJECT_SOURCE_DIR}/waveform/waveformAverage.cc ${sources} ${headers})
target_link_libraries(waveformAverage ${ROOT_LIBRARIES} Vector)

#IN FOLDER GAIN ################################################################
add_executable(makeHistDirect ${PROJECT_SOURCE_DIR}/gain/makeHistDirect.cc ${sources} ${headers})
target_link_libraries(makeHistDirect ${ROOT_LIBRARIES} Vector CustomClasses)

add_executable(makeHistIndirect ${PROJECT_SOURCE_DIR}/gain/makeHistIndirect.cc ${sources} ${headers})
target_link_libraries(makeHistIndirect ${ROOT_LIBRARIES} Vector CustomClasses)

add_executable(fitDirectLight ${PROJECT_SOURCE_DIR}/gain/fitDirectLight.cc ${sources} ${headers})
target_link_libraries(fitDirectLight ${ROOT_LIBRARIES} Vector CustomClasses)

#IN FOLDER DARKRATE ############################################################
add_executable(makePulseHist ${PROJECT_SOURCE_DIR}/darkrate/makePulseHist.cc ${sources} ${headers})
target_link_libraries(makePulseHist ${ROOT_LIBRARIES} Vector)

# Send some variables to path
set(ENV{PATH} "${PROJECT_SOURCE_DIR}/build/$ENV{PATH}")
set(ENV{LD_LIBRARY_PATH} "${PROJECT_SOURCE_DIR}/build/$ENV{LD_LIBRARY_PATH}")
