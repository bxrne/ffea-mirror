# Run only the tests I have written
import os

ffea_build_dir = os.environ["FFEA_BUILD"]
tests = ["line_connecting_rod_elements",
         "test_rngStream"]

os.chdir(ffea_build_dir)
for t in tests:
    os.system("ctest -V -R " + str(t))
    print("")
