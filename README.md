# efg-mems
MEM finding on (elastic founder) graphs

## installation

This tool uses sdsl-lite, br-index, and BDBWT.

```
git submodule update --init --recursive
cd sdsl-lite
./install.sh .
cd ..
cmake .
make
./efg-mems
```
The last command gives instructions how to use it.
You can try out the example files and shell scripts below as instructed below.
```
cd ..
wget www.cs.helsinki.fi/group/gsa/efg-mems/covid19-ecoli-efg.zip
unzip covid19-ecoli-efg.zip
cd inputs
./index-covid19.sh
./find-mems-covid19-efg.sh
./find-mems-ecoli.sh
```
To compare the results to MEM finding on br-index, you can continue as follows by installing br-index-mems in the directory containing this repository:
```
cd ..
git clone --recursive https://github.com/algbio/br-index-mems.git
cd br-index-mems
mkdir build
cd build
cmake ..
make
cd ..
cd ..
cd inputs
./find-mems-covid19-text.sh
```

## TODO
 - add br-index-mems as git submodule
 - fix bdbwt crashes/output
