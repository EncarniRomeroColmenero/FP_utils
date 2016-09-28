## Installation

These installation instructions assume you have an Ubuntu or Mac OS with the usual applications (such as ssh and a browser), but not necessarily Python or IRAF.

### Installation on Mac OS

Download the Anaconda2 installer from [Continuum Analyics' download page](https://docs.continuum.io/anaconda/install). Double-click on the pkg file and go through the installation process. Anaconda will be installed in a folder `anaconda` in your home directory.

To test whether the installation has been successful, open a new terminal window and run

```bash
python --version
```

The output should mention Anaconda, and the version should be a Python 2 (not Python 3) one.

The easiest way to install IRAF and PyRAF is to use AstroConda, which is maintained by the Space Telescope Science Institute. Add the respective channel to your Anaconda configuration.

```bash
conda config --add channels http://ssb.stsci.edu/astroconda
```

Create a new environment with the iraf, pyraf and stsci packages.

```bash
conda create -n fputil iraf pyraf stsci
```

This might take a while, so feel free to grab a cup of coffee. Once the installation is complete, activate the new environment.

```bash
source activate fputil
```

Install the pyfits library.

```bash
pip install pyfits
```

Chances are that the PyQt4 installation is broken. To check this, launch python and run

```python
import PyQt4 import QtGui
```

If this gives you an error about a missing dylib file, try the following.

```bash
cd ~/anaconda/envs/fputil/lib/python2.7/site-packages/PyQt4
install_name_tool -add_rpath /Users/christian/anaconda/lib QtGui.so
```

Clone the PySALT repository.

```bash
cd /some/convient/directory
git clone https://github.com/saltastro/pysalt
```

Create a directory `iraf` in your home directory and run the `mkiraf` command in it. Choose `xgterm` as the terminal type when prompted.

```
cd ~
mkiraf
```

A new file `login.cl` has been created in your iraf directory. Open this file in a text editor and add the following lines.

```
# Use PySALT
reset pysalt = "/some/convient/directory/pysalt/"
task pysalt.pkg = "pysalt$pysalt.cl"
reset helpdb = (envget("helpdb") // ",pysalt$lib/helpdb.mip")
```

The reset pysalt path must end with a slash. Optionally, rerun the `mkiraf` command to re-initialise the uparm parameters files, and then add the lines above to the `login.cl` file again.

You can test the PySALT installation by starting IRAF and then executing `pysalt` on the prompt. This should show you a list of the PySALT packages.

```bash
$ cl
> pysalt
 
     +----------------------------------------------------+
     |                  _ _                               |
     |   __ _ _  __ __ | | |_    http://www.salt.ac.za    |
     |  | _\ | || _| _`| |  _|   Southern African Large   |
     |  |  /__ |__|\_,_|_|\__|   Telescope PyRAF Package  |
     |  |_| \__'                                          |
     |               Development PRERELEASE               |
     |               Version 0.50  1 Oct 2014             |
     |               Recommend IRAF 2.14/PyRAF 1.8.1      |
     |               Bug reports: salthelp@salt.ac.za     |
     +----------------------------------------------------+
 
     Setting imtype=fits
 
      proptools.  saltfp.     salthrs.    saltred.    saltspec.   slottools.
```

If IRAF prompts you for curdir, you should enter the path to your IRAF directory, including a trailing slash.

You can now clone the repository for the FP utilities.

```bash
cd /some/convenient/directory
git clone git@github.com:EncarniRomeroColmenero/FP_utils.git
```

In case development should happen in a branch other than master in the remote repository, you have to create a local branch corresponding to it.

## Installation on Ubuntu

AstroConda, which is used for the instructions below, is only available for 64 bit systems.

Fist install Anaconda2. There exists no graphical Anaconda installer for Linux, but you may download a command line installer from [Continuum Analyics' download page](https://docs.continuum.io/anaconda/install).

The installer comes as a 32 and a 64 bit version. To figure out which one you need, run

```bash
uname -m
```

If the output is something like "x86_64" you need the 64 bit version, if it is something like "i686" you need the 32 bit version.

Run the downloaded installer as follows,

```bash
bash /path/to/download/anaconda-installer.sh
```

where `anaconda-installer.sh` must be replaced with the actual name of the downloaded file, of course. When prompted, choose to have Anaconda installed in your home directory and opt to prepend the Anaconda2 install location to PATH in your .bashrc file.

The rest of the installation proceeds as in the case of Mac OS, as described in the previous section. (But note that the name of the Anaconda directory might differ.)

A potential pitfall might be that a library for Qt4 is missing, in which case running the import statement

```python
from PyQt4 import QtGui
```

in Python gives an error of the form

```
ImportError: libQtGui.so.4: cannot open shared object file: No such file or directory
```

In this case install the missing library by executing

```bash
sudo apt-get install libqtgui4
```

## Running scripts on the command line

When using the command line, you have to activate the Conda environment before running any script.

```bash
source activate fputil
```

If you haven't added Anaconda to your PATH variable, use

```bash
source ~/anaconda/bin/activate fputil
```

## Running scripts in an IDE

When running a script in a IDE, you need to tell your IDE to use `~/anaconda/envs/fputil/bin/python` as the Python IDE.

While activating the fputil environment, Anaconda runs some shell scripts which set environment variables. As the IDE won't perform this activation step, you have to define these variables yourself. The easiest is to add the following code to your Bash configuration file.

```bash
# fputils Conda environment
export CONDA_PREFIX=$HOME/anaconda/envs/fputil
export CONDA_ENV_PATH=$CONDA_PREFIX
source $CONDA_ENV_PATH/etc/conda/activate.d/webbpsf-data.sh
source $CONDA_ENV_PATH/etc/conda/activate.d/iraf.sh
```







