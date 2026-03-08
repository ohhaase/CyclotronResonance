# Cyclotron Resonance Simulation Code

To build, we first need to get CMake initialized. Run

```
cmake -G "Unix Makefiles" -B build -DCMAKE_CXX_COMPILER='g++'
```

to make a build directory, force cmake to make a Unix Makefile (necessary for my windows dev environment, because it tries to compile for Visual Studio each time), and set the compiler to GNU++. 

Note: The GNU++ specification is probably no longer necessary, since we specify OpenMP with CMake. Before it was used since Clang (on Mac dev environment) didn't have OpenMP. 

Then, to generate the executable, use

```
cmake --build build
```

This generates a CyclotronResonance executable in ./build/.

To run the executable, use

```
./build/CyclotronResonance.exe <Nparticles> <Nbins> <simType> <Nthreads>
```

This will generate the resulting data in the main directory (TODO: Change output location).