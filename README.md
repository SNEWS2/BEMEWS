# BEMEWS: Better Earth Matter Effect calculation With SNEWPY

![image](https://github.com/SNEWS2/BEMEWS/assets/44247426/22dc06cb-e3b9-49d8-90ce-807a3308930f)
Image by dreamstudio.ai

## Description

`BEMEWS` is a python module for calculating the Earth-matter effect on neutrino flavor transformations. It is a standalone module that [SNEWPY](https://github.com/SNEWS2/snewpy) uses to compute the Earth-Matter Effect for supernova neutrinos. The `BEMEWS_example.py` script shows how to use the module in standalone mode. The `EarthMatter` flavor transformation class in [SNEWPY](https://github.com/SNEWS2/snewpy) is essentially the same script but with the output options turned off. 

## Installation

`BEMEWS` is available as a PyPI package and can be installed using
```
pip install BEMEWS
```
As of mid-August 2024, installation with `pip` is only available for Mac users.

If you need or wish to manually install `BEMEWS`, you can follow the instructions below:

1) You will need the packages python-devel (in Linux), [pybind11](https://pypi.org/project/pybind11/), and [setuptools](https://pypi.org/project/setuptools/)

2) Modify `setup.py` to use the correct libraries and paths. 

3) To compile enter 

`sudo python3 setup.py install`

4) If you don't want to sudo you may want to use the option

`--install-lib=destination/directory/`

## Using BEMEWS

1) The python code `BEMEWS.py` uses the module to compute the Earth matter effects upon a neutrino signal from Betelgeuse in SuperK. It will generate a lot of files in `out` folder named in the script. The output can be switched off by changing the `outputflag` parameter to `False`.

2) A script is provided that allows [SNEWPY](https://github.com/SNEWS2/snewpy) to include the Earth-matter Effect in its flavor transformation prescription. 

## Troubleshooting

1) When using `BEMEWS` you may have to set the `PYTHONPATH` environment variable to your PWD and/or wherever the `BEMEWS` module was installed in steps 3) or step 4)

2) If your script still cannot find the module you may need to put the *.so library in the same directory as the file. The *.so library is in one of the subfolders in the build directory. 

3) BEMEWS uses [OpenMP](https://www.openmp.org/). You may want to set the `OMP_NUM_THREADS` environment variable to a number suitable for your machine.

4) During compilation, you may need to set the `LIBOMP_INCLUDE` variable to enable the build. If you are a homebrew user on Mac OS X, you can do this with the command
```
export LIBOMP_INCLUDE=`brew --prefix libomp`/include
```


