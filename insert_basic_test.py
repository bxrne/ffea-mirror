import os

test_dir = input("Enter CMake test directory in allcaps (ALLCAPS): ")
test_name = input("Enter name of unit test (snake_case): ")
prefix = "/home/ryan/Software/ffea/"
path_to_test = f"tests/rods/unit/{test_name}"
cmake_text = """# 
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

set ({0} "${{PROJECT_BINARY_DIR}}/{1}/")
file (COPY {2}.ffeatest DESTINATION ${{{0}}})
add_test(NAME {2} COMMAND ${{PROJECT_BINARY_DIR}}/src/ffea {2}.ffeatest)
""".format(test_dir, path_to_test, test_name)

go = input(f"Unit test template will be written in {prefix + path_to_test}, proceed? (y/n)")
if go != "y" and go != "Y":
    quit()

try:
    os.mkdir(prefix + path_to_test)
except FileExistsError:
    print("Directory exists. Proceeding")
finally:
    os.chdir(prefix + path_to_test)

print(f"Writing {prefix}{path_to_test}/CMakeLists.txt")
with open(f"{prefix}{path_to_test}/CMakeLists.txt", "w") as f:
    f.write(cmake_text)
    
print(f"Writing {prefix}{path_to_test}/{test_name}.ffeatest")
with open(f"{prefix}{path_to_test}/{test_name}.ffeatest", "w") as f:
    f.write(test_name)
    
print(f"Appending to {prefix}{path_to_test}/../CMakeLists.txt")
with open(f"{prefix}{path_to_test}/../CMakeLists.txt", "a") as f:
    f.write(f"add_subdirectory({test_name})\n")
    
print("Done")
