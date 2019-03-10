#project name 
PROJECT(tpmsGen)

#requirement of cmake vesion
cmake_minimum_required(VERSION 3.0.0)

# Suppress MSVC CRT warnings.
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    #add_definitions(/Za)
    add_definitions(/MP) # multi-core builds
endif()

include(cmake/DownloadProject.cmake)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#Add include
set( INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include )
file(GLOB INCLUDE_DIRS ${INCLUDE_DIR}/*.h ${INCLUDE_DIR}/*.cpp )
add_library(tpmsgen ${INCLUDE_DIRS})
 target_include_directories（tpmsgen PUBLIC ${INCLUDE_DIR}）

# Eigen
if(TARGET Eigen3::Eigen)
  # If an imported target already exists, use it
  target_link_libraries(tpmsgen INTERFACE Eigen3::Eigen)
else()
	download_project(PROJ    eigen
		 			 SOURCE_DIR   ${CMAKE_CURRENT_LIST_DIR}/external/eigen
					 DOWNLOAD_DIR ${CMAKE_CURRENT_LIST_DIR}/external/.cache/eigen 
					 URL           http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz
					 URL_MD5       8ad10ac703a78143a4062c9bda9d8fd3
					
		)
  target_include_directories(tpmsgen PUBLIC ${CMAKE_CURRENT_LIST_DIR}/external/eigen)
endif()

 #Add your project files
 file(GLOB HEADFILES *.h)
file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${HEADFILES} ${SRCFILES} )
target_include_directories(${PROJECT_NAME}_bin PUBLIC ${INCLUDE_DIR} )
target_link_libraries(${PROJECT_NAME}_bin  tpmsgen)
