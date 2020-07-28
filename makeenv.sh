# GoogleTestの導入
version=1.10.0
wget https://github.com/google/googletest/archive/release-${version}.tar.gz
tar -zxvf release-${version}.tar.gz
cd googletest-release-${version}/googletest/
mkdir build
cd build
cmake ..
make
cd ../../../
rm release-${version}.tar.gz
