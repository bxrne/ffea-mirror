#!/bin/bash
# Author: Ryan Cocking (2022)
# This script builds FFEA

if [ $1 == "-h" ] || [ $1 == "--help" ];
then
    echo ""
    echo "Arguments (0|1, default 0)"
    echo "--------------------------"
    echo "do_cmake:       Build from scratch"
    echo "do_install:     Build executable after making"
    echo "do_ctest:       Enable unit testing"
    echo "restrict_tests: Only perform specific tests"
    echo "debug:          Print a lot of extra information"
    echo ""
    echo "Leave all as zero if you only want to make"
    echo ""
    exit
fi

# ============= Parameters ============= #
src_dir=$FFEA_SRC
cmake_install_prefix=$FFEA_BUILD  # not sure if needed when we're already moving to that directory
python_exec=$CONDA_PREFIX/bin/python  # if you're in the correct (python2.7) conda env, this may not be needed either
boost_root_dir=$SOFTWARE_HOME/boost_build  # leave empty for FFEA to auto-download
eigen3_include_dir=$SOFTWARE_HOME/eigen-3.3.9/include/eigen3  #  "   "   "    "

# ================ Args ================ #
do_cmake=${1:-0}  # required if making changes to the testing code
do_install=${2:-0}
do_ctest=${3:-0}
restrict_tests=${4:-0}
debug=${5:-0}
if [ $debug -eq 1 ];
then
    use_ffea_debug="ON"
else
    use_ffea_debug="OFF"
fi

echo ""
echo "do_cmake        "$do_cmake
echo "do_install      "$do_install
echo "do_ctest        "$do_ctest
echo "restrict_tests  "$restrict_tests
echo "debug           "$debug
echo ""
echo "SRC                   "$src_dir
echo "CMAKE_INSTALL_PREFIX  "$cmake_install_prefix
echo "PYTHON_EXECUTABLE     "$python_exec
echo "BOOST_ROOT            "$boost_root_dir
echo "EIGEN3_INCLUDE_DIR    "$eigen3_include_dir
echo "USE_DEBUG             "$use_ffea_debug

if [[ -z $boost_root_dir ]];
then
    echo "Building with internal Boost"
fi
if [[ -z $eigen3_include_dir ]];
then
    echo "Building with internal Eigen3"
fi                                      
# ============== Threads =============== #
max_proc_id=$(cat /proc/cpuinfo | grep 'processor' | tail -1 | awk '{print $3}')
max_threads=$((max_proc_id+1))
use_threads=$((3*max_threads/4))

# Set use_threads to 1 if empty
if [[ -z $use_threads ]];
then
    use_threads=1
fi
echo "Building with "$use_threads" threads out of "$max_threads
echo ""
# =============== Script =============== #
if [ $debug -eq 1 ];
then
    # edit the rod debug flag in-place to be switched on
    sed -i 's/dbg_print = false/dbg_print = true/g' $FFEA_SRC/src/rod_math_v9.cpp
    echo "Rod debug flag enabled using sed"
else
    sed -i 's/dbg_print = true/dbg_print = false/g' $FFEA_SRC/src/rod_math_v9.cpp
    echo "Rod debug flag disabled using sed"
fi

if [[ $do_cmake -eq 1 ]];
then
    rm -r $cmake_install_prefix
    mkdir $cmake_install_prefix
    cd $cmake_install_prefix

    if [[ -z $boost_root_dir && -z $eigen3_include_dir ]];
    then
        cmake $src_dir -DCMAKE_INSTALL_PREFIX=$cmake_install_prefix -DPYTHON_EXECUTABLE=$python_exec -DUSE_EIGEN3_INTERNAL=ON -DUSE_BOOST_INTERNAL=ON
    elif [[ -z $boost_root_dir && ! -z $eigen3_include_dir ]];
    then
        cmake $src_dir -DCMAKE_INSTALL_PREFIX=$cmake_install_prefix -DPYTHON_EXECUTABLE=$python_exec -DUSE_EIGEN3_INTERNAL=OFF -DUSE_BOOST_INTERNAL=ON -DEIGEN3_INCLUDE_DIR=$eigen3_include_dir
    elif [[ ! -z $boost_root_dir && -z $eigen3_include_dir ]];
    then
        cmake $src_dir -DCMAKE_INSTALL_PREFIX=$cmake_install_prefix -DPYTHON_EXECUTABLE=$python_exec -DUSE_EIGEN3_INTERNAL=ON -DUSE_BOOST_INTERNAL=OFF -DBOOST_ROOT=$boost_root_dir
    else
        cmake $src_dir -DCMAKE_INSTALL_PREFIX=$cmake_install_prefix -DPYTHON_EXECUTABLE=$python_exec -DUSE_EIGEN3_INTERNAL=OFF -DUSE_BOOST_INTERNAL=OFF -DBOOST_ROOT=$boost_root_dir -DEIGEN3_INCLUDE_DIR=$eigen3_include_dir -DUSE_DEBUG=$use_ffea_debug
    fi

    make -j $use_threads
else
    cd $cmake_install_prefix
fi
	
if [[ $do_install -eq 1 ]];
then
	make install -j $use_threads
else
    make -j $use_threads
fi

if [[ $do_ctest -eq 1 ]];
then
	if [[ $restrict_tests -eq 1 ]];
	then
		ctest -j $use_threads --verbose -R rod_neighbour_list_construction
	else
		ctest -j $use_threads
	fi
fi
echo "Done"
