![image](https://github.com/SNEWS2/BEMEWS/assets/44247426/22dc06cb-e3b9-49d8-90ce-807a3308930f)
Image by dreamstudio.ai


# BEMEWS
A python module for calculating the Earth-matter Effect. Works standalone or with SNEWPY

1) You will need the python-devel, pybind11 and setuptools packages

2) Modify setup.py to use the correct libraries and paths. 

3) To compile enter 

sudo python3 setup.py install 

4) If you don't want to sudo you may want to use the option

--install-lib=destination/directory/

5) The python code BEMEWS.py uses the module to compute the Earth matter effects upon a neutrino signal
   from Betelgeuse in SuperK. It will generate a lot of files in the whichever folder is named
   in the script. The output can be switched off by changing the outputflag to False

6) A script is provided that allows SNEWPY to include the Earth-matter Effect in its flavor transformation
   prescription. 

TROUBLESHOOTING:

1) When using BEMEWS you may have to set the PYTHONPATH environment variable to your PWD
   and/or wherever the BEMEWS module was installed in steps 3) or step 4)

2) If your script still cannot find the module you may need to put the *.so library in the same directory
   as the file. The *.so library is in one of the subfolders in the build directory. 

3) BEMEWS uses OpenMP. You may want to set the OMP_NUM_THREADS environment variable to a number suitable for your machine.


