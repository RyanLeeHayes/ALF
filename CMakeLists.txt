cmake_minimum_required(VERSION 3.20)
project(ALF_toplevel)

add_subdirectory(alf/bin)

# ---------------------------------------------------------------------------
# Install native binaries to <prefix>/libexec/alf/
# ---------------------------------------------------------------------------
install(TARGETS
  linear2018 nonlinear2024 whamweight impcons
  Filter Moment LM PLM
  GetLambda GetSteps
  RUNTIME DESTINATION libexec/alf
  LIBRARY DESTINATION libexec/alf
)

# ---------------------------------------------------------------------------
# Stage the Python package into the build tree with hardcoded binary paths
# ---------------------------------------------------------------------------
set(ALF_STAGE_DIR "${CMAKE_BINARY_DIR}/python")
set(ALF_LIBEXEC_DIR "${CMAKE_INSTALL_PREFIX}/libexec/alf")

file(COPY
  "${CMAKE_SOURCE_DIR}/alf"
  "${CMAKE_SOURCE_DIR}/pyproject.toml"
  "${CMAKE_SOURCE_DIR}/README.md"
  DESTINATION "${ALF_STAGE_DIR}"
)

# Drop the bin/ subtree — binaries are installed separately
file(REMOVE_RECURSE "${ALF_STAGE_DIR}/alf/bin")

# Patch Python files so they reference the installed binary location
# instead of a __file__-relative bin/ directory.
# Pattern B (longer) must be replaced before Pattern A to avoid partial match.
set(_PATCH_FILES
  loss.py
  GetLambda.py
  GetSteps.py
  postprocessDCA.py
)

foreach(_f IN LISTS _PATCH_FILES)
  set(_path "${ALF_STAGE_DIR}/alf/${_f}")
  file(READ "${_path}" _content)
  string(REPLACE
    "os.path.dirname(os.path.abspath(__file__))+'/bin/"
    "'${ALF_LIBEXEC_DIR}/"
    _content "${_content}")
  string(REPLACE
    "os.path.dirname(__file__)+'/bin/"
    "'${ALF_LIBEXEC_DIR}/"
    _content "${_content}")
  file(WRITE "${_path}" "${_content}")
endforeach()

# ---------------------------------------------------------------------------
# At install time, pip-install the staged package into the active environment
# ---------------------------------------------------------------------------
install(CODE "
  execute_process(
    COMMAND pip install \"${ALF_STAGE_DIR}\"
    RESULT_VARIABLE _pip_rc
  )
  if(NOT _pip_rc EQUAL 0)
    message(FATAL_ERROR \"pip install failed (exit code: \${_pip_rc})\")
  endif()
")
