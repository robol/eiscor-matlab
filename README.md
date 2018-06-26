# MATLAB interface to EISCOR

This repository contains a simple MATLAB interface to the polynomial
rootfinder based on unitary-plus-rank-1 QR / QZ code available at 
[https://github.com/eiscor/eiscor/](https://github.com/eiscor/eiscor). 

To compile this mex file you need to install a copy of eiscor before. 
You can download and install eiscor with the following commands: 

```Bash
git clone git://github.com/eiscor/eiscor
cd eiscor && make && make install
```

If the compilation and installation succeeds, then you can compile this
mex file running the following command at the MATLAB prompt. 

```Matlab 
mex eiscor_roots.F90 ~/eiscor/lib/libeiscor.so.0.2.0
````

   
