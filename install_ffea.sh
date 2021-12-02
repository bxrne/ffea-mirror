#!/bin/bash
doinstall=${2:0}
dotest=${3:0}
if [[ $doinstall -eq 1 ]];
then
	rm -r $FFEA_BUILD
	mkdir $FFEA_BUILD
	cd $FFEA_BUILD
	cmake ../ffea -DCMAKE_INSTALL_PREFIX=$FFEA_BUILD -DPYTHON_EXECUTABLE=/home/ryan/Software/anaconda3/envs/py2/bin/python
	make -j2
	make install -j2
elif [[ $doinstall -eq 0 ]];
then
	cd $FFEA_BUILD
	make -j2
fi

if [[ $dotest -eq 1 ]];
then
	ctest -j2
fi
