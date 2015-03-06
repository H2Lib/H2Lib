find_package(Doxygen)
if(DOXYGEN_FOUND)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/Doc)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
else()
  message(STATUS "Doxygen not found -- not building docs")
endif()

if(NOT WIN32)
  add_custom_target(cleandoc
      COMMAND rm -rf ${PROJECT_SOURCE_DIR}/Doc/html ${PROJECT_SOURCE_DIR}/Doc/latex
      COMMENT "Deleting documentation."
      )
endif()

