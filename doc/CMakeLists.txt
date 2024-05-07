file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html)

file(GLOB_RECURSE allDocFiles ${CMAKE_SOURCE_DIR}/*.cpp *.h)
string(REGEX REPLACE "${CMAKE_BINARY_DIR}[^;]+;?" "" allDocFiles "${allDocFiles}")

set(OUTPUT_DOC_DIR ${CMAKE_CURRENT_BINARY_DIR}/)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.cmake ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

set(targetName "doc_${PROJECT_NAME}")

add_custom_command(

    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doc.txt
    OUTPUT ${targetName}

    COMMENT "doxygen documentation generation for ${PROJECT_NAME}"

    COMMAND doxygen ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    DEPENDS ${allDocFiles}

    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set(OthertargetName "mkdocs_${PROJECT_NAME}")

# ######################################################
# ######################################################
# ######################################################
# Set the docs directory path
set(DOCS_DIR "${CMAKE_CURRENT_BINARY_DIR}/docs/")
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/docs/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/docs/)

# Build the list of all *.md files found in tests directory
set(EXAMPLES_DIR "${CMAKE_SOURCE_DIR}/tests")
file(GLOB_RECURSE EXAMPLE_FILES "${EXAMPLES_DIR}/description.md")
set(EXAMPLE_LIST "")

set(DOCS_EXAMPLES_DIR "Started/Examples")

foreach(EXAMPLE_FILE ${EXAMPLE_FILES})
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_FILE} NAME_WE)
    get_filename_component(EXAMPLE_FILE_DIR ${EXAMPLE_FILE} DIRECTORY)
    string(REPLACE "${EXAMPLES_DIR}" "" EXAMPLE_FILE_DIR ${EXAMPLE_FILE_DIR})
    string(REPLACE "tests" "" EXAMPLE_FILE_DIR ${EXAMPLE_FILE_DIR})
    string(REPLACE "/" " " EXAMPLE_NAV_NAME ${EXAMPLE_FILE_DIR})
    configure_file(${EXAMPLE_FILE} "${DOCS_DIR}/${DOCS_EXAMPLES_DIR}/${EXAMPLE_FILE_DIR}/${EXAMPLE_NAME}.md" COPYONLY)
    set(EXAMPLE_LIST "${EXAMPLE_LIST}\n      - ${EXAMPLE_NAV_NAME}: ${DOCS_EXAMPLES_DIR}/${EXAMPLE_FILE_DIR}/${EXAMPLE_NAME}.md")
endforeach()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mkdocs.yml.cmake ${CMAKE_CURRENT_BINARY_DIR}/mkdocs.yml @ONLY)
add_custom_command(

    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/mkdocs.txt
    OUTPUT ${OthertargetName}

    COMMENT "mkdocs documentation generation for ${PROJECT_NAME}"

    COMMAND mkdocs build
    DEPENDS ${mkConfile}

    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(${targetName} DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doc.txt)
add_custom_target(${OthertargetName} DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/mkdocs.txt ${targetName})
add_dependencies(doc ${OthertargetName})


# install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/site DESTINATION doc_${PROJECT_NAME})
