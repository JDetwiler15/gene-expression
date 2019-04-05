v=3.7.2
py_url=https://www.python.org/ftp/python/$v/Python-$v.tar.xz
tar_name=$(basename $py_url)
dir_name=~/local/Python-$v

echo "Downloading $tar_name from python.org..."
wget $py_url 2>/dev/null
echo "Untarring $tar_name..."
tar Jxf $tar_name
echo "Relocating it to $dir_name"
mkdir -p $dir_name
rm -f $tar_name
mv Python-$v ~/local
cd $dir_name

echo "Now configuring it to be installed locally."
./configure --prefix=$dir_name 1>/dev/null
echo "Compiling Python-$v with -j 40."
make -j 40 &>/dev/null
echo "Building binaries with -j 40."
make install -j 40 &>/dev/null

echo "Linking Python and Pip binaries to ~/bin."
mkdir -p ~/bin
cd ~/bin
ln -s $dir_name/bin/python${v:0:1} python${v:0:1}
ln -s $dir_name/bin/pip${v:0:1} pip${v:0:1}
echo "You will want to prepend PATH with ~/bin to use this Python install."
echo "Edit your ~/.bashrc with 'export PATH=~/bin:\$PATH'."
echo "Then, run install-packages.sh to grab NetworkX and others."
