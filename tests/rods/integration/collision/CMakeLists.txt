#
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

set (COLLISION_DIR          "${PROJECT_BINARY_DIR}/tests/rods/integration/collision/")
set (COLLISION_TEST_NAME    int_rodStericCollision_pushApart)

file (COPY run_test.py            DESTINATION ${COLLISION_DIR})
file (COPY ${COLLISION_TEST_NAME}.ffeatest  DESTINATION ${COLLISION_DIR})

# add_test(NAME ${COLLISION_TEST_NAME} COMMAND ${PROJECT_BINARY_DIR}/src/ffea ${COLLISION_TEST_NAME}.ffeatest)
add_test(NAME ${COLLISION_TEST_NAME} COMMAND ${PYTHON_EXECUTABLE} run_test.py)
# set_tests_properties(${COLLISION_TEST_NAME} PROPERTIES ENVIRONMENT PYTHONPATH=${PROJECT_SOURCE_DIR}/ffeatools/modules:$ENV{PYTHONPATH})
