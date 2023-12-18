# efg-mems
MEM finding on (elastic founder) graphs

## installation

This tool requires sdsl-lite, br-index, and BDBWT.
Latter two should be installed in the parent directory of efg-mems.

```
mkdir efg-mems-root
cd efg-mems-root
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh
cd ..
git clone --recursive https://github.com/U-Ar/br-index.git
git clone --recursive https://github.com/algbio/bdbwt
git clone https://github.com/algbio/efg-mems
cd efg-mems
cmake .
make
./efg-mems
```
The last command gives instructions how to use it.
You can try out the example files and shell scripts below as instructed below.
For these to work you need to install br-index (check a closed issue for tips).
```
cd ..
wget www.cs.helsinki.fi/group/gsa/efg-mems/covid19-ecoli-efg.zip
unzip covid19-ecoli-efg.zip
cd inputs
./index-covid19.sh
./find-mems-covid19-efg.sh
./find-mems-ecoli.sh
```
To compare the results to MEM finding on br-index, you can continue as follows:
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
 - parse input queries with the correct FASTA header info
 - finalize linear MEMS file format
 - output both formats (GAF/MEMS) in one execution
