# 

set (STERICENERGYDIR "${PROJECT_BINARY_DIR}/tests/rods/unit/steric_energy_two_rod_elements/")
file (COPY steric_energy_two_rod_elements.ffeatest DESTINATION ${STERICENERGYDIR})
add_test(NAME steric_energy_two_rod_elements COMMAND ${PROJECT_BINARY_DIR}/src/ffea steric_energy_two_rod_elements.ffeatest)
