# efg-mems
MEM finding on (elastic founder) graphs

# installation

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
.\efg-mems
```
The last command gives instructions how to use it.
You can try out the example files and shell scripts:
```
cd ..
wget www.cs.helsinki.fi/u/vmakinen/ef-mems/covid19-ecoli-efg.zip
unzip covid19-ecoli-efg.zip
cd inputs
./find-mems-covid19.sh
./find-mems-covid19-2.sh
./find-mems-ecoli.sh
```
