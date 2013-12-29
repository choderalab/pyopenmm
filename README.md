
pycuda auto-install issues:
```
Installed /Users/choderaj/code/pyopenmm-experimental/install/lib/python2.7/site-packages/pyopencl-2013.2-py2.7-macosx-10.6-x86_64.egg
Searching for pycuda
Reading https://pypi.python.org/simple/pycuda/
Best match: pycuda 2013.1.1
Downloading https://pypi.python.org/packages/source/p/pycuda/pycuda-2013.1.1.tar.gz#md5=acf9319ab2970d9700ed6486aa87b708
Processing pycuda-2013.1.1.tar.gz
Writing /var/folders/m1/wl3ptzqx2csf2r_q2bfy78n0h7kcfv/T/easy_install-GYYjwW/pycuda-2013.1.1/setup.cfg
Running pycuda-2013.1.1/setup.py -q bdist_egg --dist-dir /var/folders/m1/wl3ptzqx2csf2r_q2bfy78n0h7kcfv/T/easy_install-GYYjwW/pycuda-2013.1.1/egg-dist-tmp-gJOhhb
*************************************************************
*** I have detected that you have not run configure.py.
*************************************************************
*** Additionally, no global config files were found.
*** I will go ahead with the default configuration.
*** In all likelihood, this will not work out.
*** 
*** See README_SETUP.txt for more information.
*** 
*** If the build does fail, just re-run configure.py with the
*** correct arguments, and then retry. Good luck!
*************************************************************
*** HIT Ctrl-C NOW IF THIS IS NOT WHAT YOU WANT
*************************************************************
Continuing in 1 seconds...    
warning: no files found matching '*.cpp' under directory 'bpl-subset/bpl_subset/boost'
warning: no files found matching '*.html' under directory 'bpl-subset/bpl_subset/boost'
warning: no files found matching '*.inl' under directory 'bpl-subset/bpl_subset/boost'
warning: no files found matching '*.txt' under directory 'bpl-subset/bpl_subset/boost'
warning: no files found matching '*.h' under directory 'bpl-subset/bpl_subset/libs'
warning: no files found matching '*.ipp' under directory 'bpl-subset/bpl_subset/libs'
warning: no files found matching '*.pl' under directory 'bpl-subset/bpl_subset/libs'
Compiling with an SDK that doesn't seem to exist: /Developer/SDKs/MacOSX10.6.sdk
Please check your Xcode installation
clang: warning: no such sysroot directory: '/Developer/SDKs/MacOSX10.6.sdk'
In file included from src/cpp/cuda.cpp:1:
In file included from src/cpp/cuda.hpp:12:
/usr/local/cuda/include/cuda.h:53:10: fatal error: 'stdlib.h' file not found
#include <stdlib.h>
         ^
1 error generated.
error: Setup script exited with error: command 'gcc' failed with exit status 1
```