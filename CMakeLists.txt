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
set( TPMSGEN_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include )
file(GLOB INCLUDE_DIRS ${TPMSGEN_INCLUDE_DIR}/*.h ${TPMSGEN_INCLUDE_DIR}/*.cpp )
add_library(tpmsgen ${INCLUDE_DIRS})
target_include_directories(tpmsgen PUBLIC ${TPMSGEN_INCLUDE_DIR})

# Eigen
if(TARGET Eigen3::Eigen)
  # If an imported target already exists, use it
  target_link_libraries(tpmsgen INTERFACE Eigen3::Eigen)
else()
	download_project(PROJ    eigen
		 			 SOURCE_DIR   ${CMAKE_CURRENT_LIST_DIR}/external/eigen
					 DOWNLOAD_DIR ${CMAKE_CURRENT_LIST_DIR}/external/.cache/eigen 
<<<<<<< HEAD
					 GIT_REPOSITORY http://bitbucket.org/eigen/eigen/get/3.7.7.tar.gz
					 GIT_TAG 3.3.7
=======
					 URL           http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
					 #URL_MD5       8ad10ac703a78143a4062c9bda9d8fd3
>>>>>>> 5484bc40cf5bddf18dfaed605910d67411e8540d
					
		)
  target_include_directories(tpmsgen PUBLIC ${CMAKE_CURRENT_LIST_DIR}/external/eigen)
endif()

 #Add your project files
 file(GLOB HEADFILES *.h)
file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${HEADFILES} ${SRCFILES} )
target_include_directories(${PROJECT_NAME}_bin PUBLIC ${TPMSGEN_INCLUDE_DIR} )
target_link_libraries(${PROJECT_NAME}_bin  tpmsgen)
