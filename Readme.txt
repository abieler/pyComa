Install Anaconda:
*****************
1)  download python anaconda distribution from http://continuum.io/downloads

1a) execute "bash <downloaded file>" in the shell.

1b) choose installation path of anaconda.


Config Anaconda:
****************
1c) execute "conda update conda" in the shell.

1d) execute "conda update anaconda" in the shell.

1e) execute "conda install mpi4py" in the shell.

1f) execute "conda install bokeh" in the shell.

1g) modify your PATH variable such that python is called
from /YourPathTo/anaconda/bin/python
(newer anaconda distros can do this automatically for you
during the install process)


Install SPICE/ pySPICE
***********************
pySPICE in an unofficial wrapper for the NAIF cspice toolkit

2a) download cspice toolkit from:
http://naif.jpl.nasa.gov/naif/toolkit_C.html.

2b) download pySpice from
https://github.com/rca/PySPICE.

2c) extract both downloaded packages and move the extracted cspice
folder into the PySPICE-master directory.

2d) from inside PySPICE-master execute the following two commands from
the shell:
python setup.py build_ext
python setup.py install


Install Git
***********


Install pyComa
**************
4a) execute "git clone https://github.com/abieler/pyComa" in the shell.
or download it manually from above url.

4b) execute "python setup.py build_ext --inplace" in the shell inside
the pyComa/bin/utils directory.

4c) create the following symlink (if needed):
sudo ln -s /YourPathTo/anaconda/bin /opt/local/anaconda/bin


Workarounds
***********
On OSX mpi4py might be broken. to fix this, create the following symlink:

sudo ln -s /yourPathTo/anaconda /opt/anaconda1anaconda2anaconda3
