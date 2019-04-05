mkdir ~/R
cd ~/R
wget https://cloud.r-project.org/src/base/R-3/R-3.5.3.tar.gz
tar -xzvf R-3.5.3.tar.gz
cd R-3.5.3
./configure --with-x=no
make
