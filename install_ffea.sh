#!/bin/bash
# Author: Ryan Cocking (2022)
# This script builds FFEA, assuming you have pre-compiled eigen3 and boost.
do_cmake=${1:-0}
do_install=${2:-0}
do_test=${3:-0}
restrict_tests=${3:-1}

src_dir=$FFEA_SRC
cmake_install_prefix=$FFEA_BUILD  # not sure if needed when we're already moving to that directory
python_exec=$CONDA_PREFIX/bin/python  # if you're in the correct (python2.7) conda env, this may not be needed either
boost_root_dir=$SOFTWARE_HOME/boost_build
eigen3_include_dir=$SOFTWARE_HOME/eigen-3.3.9/include/eigen3

echo ""
echo "do_cmake = "$do_cmake
echo "do_install = "$do_install
echo "do_ctest = "$do_ctest
echo "restrict_tests = "$restrict_tests
echo ""

# CMake is required if making changes to the testing environment
if [[ $do_cmake -eq 1 ]];
then
    rm -r $cmake_install_prefix
    mkdir $cmake_install_prefix
    cd $cmake_install_prefix
    cmake $src_dir -DCMAKE_INSTALL_PREFIX=$cmake_install_prefix -DPYTHON_EXECUTABLE=$python_exec -DUSE_EIGEN3_INTERNAL=OFF -DUSE_BOOST_INTERNAL=OFF -DBOOST_ROOT=$boost_root_dir -DEIGEN3_INCLUDE_DIR=$eigen3_include_dir
    make -j8
else
    cd $cmake_install_prefix
fi

make -j8
	
if [[ $do_install -eq 1 ]];
then
	make install -j8
fi

if [[ $do_ctest -eq 1 ]];
then
	if [[ $restrict_tests -eq 1 ]];
	then
		ctest -j8 --verbose -R rod_neighbour_list_construction
	else
		ctest -j8
	fi
fi
echo "Done"
