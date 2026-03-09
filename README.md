# Cyclotron Resonance Simulation Code

To build, we first need to get CMake initialized. Run

```
cmake -G "Unix Makefiles" -B build -DCMAKE_CXX_COMPILER='g++'
```

This makes a build directory, forces cmake to make a Unix Makefile, and sets the compiler to GNU++. 

Note: The GNU++ specification might no longer be necessary, since we specify OpenMP with CMake. Before it was used since Clang (on Mac dev environment) didn't have OpenMP. 

Then, to generate the executable, use

```
cmake --build build
```

This generates a CyclotronResonance executable in ./build/.

To run the executable, use (include `.exe` if on Windows)

```
./build/CyclotronResonance(.exe) <Nparticles> <Nbins> <simType> <Nthreads>
```

This will generate the resulting data in the main directory (TODO: Change output location).