##############################################################
###    Build Library                                       ###
##############################################################

# all local .c files are included 
FILE(GLOB srcfiles "${CMAKE_CURRENT_SOURCE_DIR}/*.c")
add_library (H2Lib SHARED ${srcfiles})
target_link_libraries( H2Lib LINK_PUBLIC ${BLAS_LIBRARIES} )

# command to specify a seperate include directory (currently all
#   include-files are stored together with the code files)
target_include_directories (H2Lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


##############################################################
###    Install Rule                                        ###
##############################################################
FILE(GLOB incfiles "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

install(TARGETS H2Lib DESTINATION lib)
install(FILES ${incfiles} DESTINATION include)

