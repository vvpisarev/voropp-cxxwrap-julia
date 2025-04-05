# voropp-cxxwrap-julia
C++ side of Julia wrapper for Voro++

## Build Instructions

### Building Voro++
1. First, build Voro++ library (https://github.com/chr1shr/voro) by cloning the repo in a local directory:
```
cd /opt
sudo clone https://github.com/chr1shr/voro.git
```
2. In voro repo folder execute the cmake build commands:

```
cd voro
mkdir build && cd build
cmake -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_INSTALL_PREFIX=../_install ..
cmake --build . --config Release
make install
```
This will install `libvoro++.a` into `$(CMAKE_INSTALL_PREFIX)/lib64`

### Configure and build this library
In the top directory of this library:
```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="/path/to/libcxxwrap-julia-prefix;/path/to/voro/_install" ../
cmake --build . --config Release
```

The `/path/to/libcxxwrap-julia-prefix` for `CMAKE_PREFIX_PATH` can be obtained from Julia using:
```julia
julia> using CxxWrap
julia> CxxWrap.prefix_path()
```

The `/path/to/voro/` is the path where you installed Voro++ on the previous step.