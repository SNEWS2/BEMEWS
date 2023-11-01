# EMEWS
A python module for calculating the Earth-matter Effect. Works standalone or with SNEWPY

1) You will need the python-devel, pybind11 and setuptools packages

2) Modify setup.EMEWS.py to use the correct libraries and paths. 

3) To compile enter 

sudo python3 setup.EMEWS.py install 

4) If you don't want to sudo you may want to use the option

--install-lib=destination/directory/

5) The file EMEWS.py uses the module to compute the Earth matter effects upon a neutrino signal
   from Betelgeuse in SuperK. It will generate a lot of files in the whichever folder is picked
   in the script. The output can be switched off by changing the outputflag to False

6) A script is provided that allows SNEWPY to include the Earth-matter Effect in its flavor transformation
   prescription. The results from EMEWS are used by SNEWPY to replace the D's that appear in the transformation formulae.  

TROUBLESHOOTING:

1) When running EMEWS.py you may have to set the PYTHONPATH environment variable to your PWD
   and/or wherever the EMEWS module was installed in steps 3) or step 4)

2) If the script still cannot find the module you may need to put the *.so library in the same directory
   as the EMEWS.py file. The *.so library is in one of the subfolders in the build directory. 

3) You may want to set the OMP_NUM_THREADS environment variable to a number suitable for your machine


