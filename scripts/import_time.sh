#! /bin/bash

echo PWD is $PWD

for path in thompson thompson/drawing thompson/examples
do
	dir="../$path/__pycache__"
	if [ -d $dir ]
	then
		rm -r $dir
	fi
done
python3 -m cProfile -s cumtime import_time.py > result.txt
head result.txt