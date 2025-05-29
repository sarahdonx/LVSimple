PROJECT(mysolver)
CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

# List your source files here - this example has just one.
SET(SOURCE_FILES ActiveContraction.C
BoundaryConditionsHV01.C    
MechanicsModelHV01.C       	     
ModelInitialization.C  			 
tetVolumeCalculation.C 		
main.C)
ADD_EXECUTABLE(main3dSimple ${SOURCE_FILES})

FIND_PACKAGE(IBAMR 0.13.0 REQUIRED)
TARGET_LINK_LIBRARIES(main3dSimple IBAMR::IBAMR3d)
TARGET_INCLUDE_DIRECTORIES(main3dSimple PUBLIC "./")
# IBAMR saves the flags it used to compile - you can reuse them if you want to
SET(CMAKE_CXX_FLAGS ${IBAMR_CXX_FLAGS})
