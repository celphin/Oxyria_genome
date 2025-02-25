#########################################
#Installing mini conda in home directory:
  #Needed in min
##########################################

curl "https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh" -o miniconda.sh
 #follow prompts
 #when asked to modify ~/bashrc  say 'no'
chmod +x miniconda.sh

~/miniconda.sh

Do you accept the license terms? [yes|no]
[no] >>>
Please answer 'yes' or 'no': yes

Miniconda2 will now be installed into this location:
/home/celphin/miniconda2

  - Press ENTER to confirm the location
  - Press CTRL-C to abort the installation
  - Or specify a different location below
#press enter


Preparing transaction: done
Executing transaction: done
installation finished.
WARNING:
    You currently have a PYTHONPATH environment variable set. This may cause
    unexpected behavior when running the Python interpreter in Miniconda2.
    For best results, please verify that your PYTHONPATH only points to
    directories of packages that are compatible with the Python interpreter
    in Miniconda2: /home/msandler/miniconda2
Do you wish the installer to initialize Miniconda2
by running conda init? [yes|no] ...
Do you wish the installer to initialize Miniconda2
by running conda init? [yes|no]
[no] >>> yes
no change     /home/celphin/miniconda2/condabin/conda
no change     /home/celphin/miniconda2/bin/conda
no change     /home/celphin/miniconda2/bin/conda-env
no change     /home/celphin/miniconda2/bin/activate
no change     /home/celphin/miniconda2/bin/deactivate
no change     /home/celphin/miniconda2/etc/profile.d/conda.sh
no change     /home/celphin/miniconda2/etc/fish/conf.d/conda.fish
no change     /home/celphin/miniconda2/shell/condabin/Conda.psm1
no change     /home/celphin/miniconda2/shell/condabin/conda-hook.ps1
no change     /home/celphin/miniconda2/lib/python2.7/site-packages/xontrib/conda.xsh
no change     /home/celphin/miniconda2/etc/profile.d/conda.csh
modified      /home/celphin/.bashrc

==> For changes to take effect, close and re-open your current shell. <==

If you'd prefer that conda's base environment not be activated on startup,
   set the auto_activate_base parameter to false:

conda config --set auto_activate_base false

Thank you for installing Miniconda2!

#-------------------------------
rm miniconda.sh
~/miniconda2/bin/conda

#------------------
usage: conda [-h] [-V] command ...

conda is a tool for managing and deploying applications, environments and packages.

Options:

positional arguments:
  command
    clean        Remove unused packages and caches.
    config       Modify configuration values in .condarc. This is modeled
                 after the git config command. Writes to the user .condarc
                 file (/home/celphin/.condarc) by default.
    create       Create a new conda environment from a list of specified
                 packages.
    help         Displays a list of available conda commands and their help
                 strings.
    info         Display information about current conda install.
    init         Initialize conda for shell interaction. [Experimental]
    install      Installs a list of packages into a specified conda
                 environment.
    list         List linked packages in a conda environment.
    package      Low-level conda package utility. (EXPERIMENTAL)
    remove       Remove a list of packages from a specified conda environment.
    uninstall    Alias for conda remove.
    run          Run an executable in a conda environment. [Experimental]
    search       Search for packages and display associated information. The
                 input is a MatchSpec, a query language for conda packages.
                 See examples below.
    update       Updates conda packages to the latest compatible version.
    upgrade      Alias for conda update.

optional arguments:
  -h, --help     Show this help message and exit.
  -V, --version  Show the conda version number and exit.

conda commands available from other packages:
  env
