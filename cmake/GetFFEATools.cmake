
set(FFEATOOLS_DOWNLOAD_VERSION "master")

include(FetchContent)
FetchContent_Declare(
  FFEATools
  GIT_REPOSITORY https://bitbucket.org/FFEA/ffeatools.git
  GIT_TAG ${FFEATOOLS_DOWNLOAD_VERSION}
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
#set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
FetchContent_MakeAvailable(FFEATools)


# FFEAtools requires Python, but doesn't have it's own CMakeLists
if(PYTHON3_EXACT_VERSION)
    set(PYTHON3_EXACT_VERSION_ARG ${PYTHON3_EXACT_VERSION} EXACT)
endif()
find_package(Python3 ${PYTHON3_EXACT_VERSION_ARG} REQUIRED COMPONENTS Interpreter) #  Development.Module not required, not making wheels?
unset(PYTHON3_EXACT_VERSION_ARG)
message(STATUS "Python found at " ${Python3_EXECUTABLE})

# Function to find if python module MODULE_NAME is available, and error if it is not available.
function(ffea_search_python_module MODULE_NAME)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import ${MODULE_NAME}; print(${MODULE_NAME}.__version__) if hasattr(${MODULE_NAME}, '__version__') else print('Unknown');"
    RESULT_VARIABLE _RESULT
    OUTPUT_VARIABLE MODULE_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  if(${_RESULT} STREQUAL "0")
    message(STATUS "Found python module: ${MODULE_NAME} (version \"${MODULE_VERSION}\")")
  else()
    message(FATAL_ERROR 
      "  Unable to find required python module \"${MODULE_NAME}\".\n"
      "  Please install this to your python environment, e.g.\n"
      "    python3 -m pip install ${MODULE_NAME}")
  endif()
endfunction()

# Create a venv and install ffeatools
# Look for python module venv, error if not found
ffea_search_python_module(venv)
set(VENV_EXECUTABLE ${Python3_EXECUTABLE} -m venv)
set(VENV_DIR ${CMAKE_BINARY_DIR}/venv)
set(VENV_BYPRODUCTS "${VENV_DIR}")
if(WIN32)
    set(VENV_PIP "${VENV_DIR}\\Scripts\\pip.exe")
else()
    set(VENV_PIP ${VENV_DIR}/bin/pip)
endif()
# make a venv to install our python package in it
add_custom_target(ffeatools_venv
# create venv
COMMAND ${VENV_EXECUTABLE} ${VENV_DIR}
# Manually install ffeatools and it's dependencies
COMMAND ${CMAKE_COMMAND} -E env CC=gcc ${VENV_PIP} install --editable "${ffeatools_SOURCE_DIR}" -v
BYPRODUCTS ${VENV_BYPRODUCTS}
WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
COMMENT "Install ${PYTHON_DISTRIBUTION_NAME} FFEATools into ${VENV_DIR}"
)