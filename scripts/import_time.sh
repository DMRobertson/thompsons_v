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
python3 -c "import sys; s1 = set(sys.modules); import thompson; s2 = set(sys.modules); import pprint; print(len(s1), len(s2)); pprint.pprint(sorted(s2 - s1))" > modules.txt
