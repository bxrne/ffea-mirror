# 

set (TWOSPHEREDIR "${PROJECT_BINARY_DIR}/tests/rods/unit/two_sphere_volume_intersection/")
file (COPY two_sphere_volume_intersection.ffeatest DESTINATION ${TWOSPHEREDIR})
add_test(NAME two_sphere_volume_intersection COMMAND ${PROJECT_BINARY_DIR}/src/ffea two_sphere_volume_intersection.ffeatest)
