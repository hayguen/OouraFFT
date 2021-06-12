
---

# General Purpose FFT (Fast Fourier/Cosine/Sine Transform) Package

---

<!-- toc -->

- [Brief Description](#brief-description)
- [Dependencies](#dependencies)
- [Git](#git)
- [CMake](#cmake)
- [Functions](#functions)
- [License](#license)

<!-- tocstop -->

---

## Brief description:

The General Purpose FFT (Fast Fourier/Cosine/Sine Transform) Package:
* ANSI C library
* **was** only double precision floating point 'double'

This is a fork of his FFT Package - currently with some changes for the 1D-FFT:
* added cmake support, to build libraries
* added support for single precision 'float'
* added include headers
* enhanced tests/samples
* use ctest

Origin of this package is [https://www.kurims.kyoto-u.ac.jp/~ooura/fft.html](https://www.kurims.kyoto-u.ac.jp/~ooura/fft.html)


## Dependencies:

On Debian/Ubuntu Linux following packages should be installed:

```
sudo apt-get install build-essential gcc cmake git
```


## Git:
This archive's source can be downloaded with git:
```
git clone https://github.com/hayguen/OouraFFT.git
```


## CMake:
There's now CMake support to build the many libraries from the source files.

There are a few CMake options to modify library size and optimization.
You can explore all available options with `cmake-gui` or `ccmake`,
the console version - after having installed (on Debian/Ubuntu Linux) one of
```
sudo apt-get install cmake-qt-gui
sudo apt-get install cmake-curses-gui
```

The options:
* `ARCH` to set the target architecture for the compiler
* `OOURA_USE_FLOAT_PREC` to compile the libraries for single precision 'float', default is 'double'
* `OOURA_USE_FAST_MATH` to (de)activate fast floating point math

Options can be passed to `cmake` at command line, e.g.
```
cmake -DARCH=core2 -DOOURA_USE_FLOAT_PREC=ON
```

My Linux distribution defaults to GCC. With installed CLANG and the bash shell, you can use it with
```
mkdir build
cd build
CC=/usr/bin/clang CXX=/usr/bin/clang++ cmake -DCMAKE_BUILD_TYPE=Release ../
cmake -DCMAKE_BUILD_TYPE=Release -DARCH=core2 -DOOURA_USE_FLOAT_PREC=ON ../
ccmake .                          # or: cmake-gui .
cmake --build .                   # or simply: make
```


## Functions:

Routines in the Package:
* cdft: Complex Discrete Fourier Transform
* rdft: Real Discrete Fourier Transform
* ddct: Discrete Cosine Transform
* ddst: Discrete Sine Transform
* dfct: Cosine Transform of RDFT (Real Symmetric DFT)
* dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)

For details on parameters, see the original [readme.txt](readme.txt)
or the new header files [include/ooura/](include/ooura/).

The 2D and 3D-routines are `as-is` from Takuya OOURA.
See his notes in [readme2d.txt](readme2d.txt)


## License:
MIT License, see [LICENSE](LICENSE)

Copyright(c) 2021 Hayati Ayguen

Copyright(C) 1996-2001 Takuya OOURA
