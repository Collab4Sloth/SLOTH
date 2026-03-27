set(Python3_FIND_STRATEGY LOCATION)
set(Python3_FIND_VIRTUALENV FIRST)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

function(check_python_module module result_var)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import ${module}"
        RESULT_VARIABLE _res
        OUTPUT_QUIET ERROR_QUIET
    )
    if(_res EQUAL 0)
        set(${result_var} TRUE PARENT_SCOPE)
    else()
        set(${result_var} FALSE PARENT_SCOPE)
    endif()
endfunction()