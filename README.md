## Installation
```bash
mamba env create -n ogtk -c conda-forge pyyaml numpy pandas matplotlib
mamba activate ogtk

git clone --recursive https://github.com/jgi-jbredeson/ogtk.git
cd ogtk
make install PREFIX=$(mamba env list | grep ortho | tr -s ' ' | cut -f3 -d' ')
```
