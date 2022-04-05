# INF443 Project

## Introduction

This project uses the CGP library, which can be found [here](https://github.com/drohmer/CGP/).

[See CGP Documentation](https://imagecomputing.net/cgp/index.html) for details.

## Compiling the project

The directory _library/_ contains the source code of CGP, while the directory _scenes/_ contains the scenes.
Each scene is an independant program with its own Makefile and/or CMakeLists.txt. 

The scenes should be run from their root path where _shaders/_ (and possibly _assets/_) directories are accessible.

### Dependencies

CGP requires
* A C++14 (or greater) compatible compiler (GCC/CLang, or a recent Visual Studio).
* An OpenGL 3.3 (or greater) compatible system.
* [libGLFW](https://www.glfw.org/) and [pkgconfig](https://www.freedesktop.org/wiki/Software/pkg-config/) installed for Linux/MacOS system.

### Linux/MacOS

Assuming a command line opened in one of the scenes.

* _Method 1._ Using the provided Makefile:
```c++
$ make
$ ./[executable-name]
```

* _Method 2._ Using the provided CMakeLists.txt:
```c++
$ mkdir build
$ cd build
$ cmake ..
$ make
$ cd ..
$ build/[executable-name]
```

### Windows


* _Method 1._ Create a Visual Studio project using CMake
* _Method 2._ Open the CMakeLists.txt using the internal CMake tool from Visual.

_Once opened by Visual Studio, the project should be configured to compile and be executed directly without further setup. Make sure your Windows version is updated for Visual Studio to be able to compile correctly C++14._


### Detailed system set-up and compilation

A detailed tutorial on how to install and compile C++ code is available here if needed: [Detailed installation and compilation for CGP](https://imagecomputing.net/cgp/compilation).
