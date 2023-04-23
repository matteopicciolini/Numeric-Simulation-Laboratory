# Laboratory of Numeric Simulation
## Matteo Picciolini
This directory contains the exercises of the Laboratory of Numeric Simulation A.A. 2022-2023.


### Prerequisites
In order to execute the comnands described in this file you will need several tools, such as `cmake`, `make`, `gcc`. On a Debian-based system you can install them with:

```
sudo apt-get install cmake make gcc
```



### Build
`cmake` introduces a configuration step, just before the build step. `cmake` input consists of [CMakeLists.txt](CMakeLists.txt) files, which describes the project to be built, `*.cmake` files, which contains collections of cmake macro that can be used in `CMakeLists.txt` for specific tasks, `*.in` files that are used as templates for automatically generated files.


### Configure and Compile ###
1. Clone the project repository onto your local machine using Git or download the zip file and extract it.

2. Open a terminal or command prompt and navigate to the root directory of the project.

3. Create a new directory called "build" (or any other name you prefer), navigate to the build directory and run CMake to generate the build files using the following commands:
```
mkdir -p build
cd build
cmake .. 
```
4. Once CMake has finished generating the build files, you can build the project using the following commands.
If you want to compile all the Exercise_XX.X.cpp files contained in all the Lezioni_XX folders, you can use:
```
make
```
Otherwise, if you want to compile a specific file, you can type:
```
cd Lesson_XX
make Exercise_XX.X.cpp
```