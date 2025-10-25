# Note that this script can accept some limited command-line arguments, run
# `julia build_tarballs.jl --help` to see a usage message.
using BinaryBuilder, Pkg

name = "voropp_wrapper"
version = v"0.1.0"

basepath = dirname(@__DIR__)
# Collection of sources required to complete build
sources = [
    GitSource("https://github.com/chr1shr/voro.git", "2cb6cefc690be1c653bfb8e65559ee8441c0b21f"),
    DirectorySource(basepath; target = "voropp-cxxwrap-julia")
]

# Bash recipe for building across all platforms
script = raw"""
cd $WORKSPACE/srcdir
cd voro
cmake -B build -DCMAKE_INSTALL_PREFIX=$prefix -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_BUILD_TYPE=Release -DVORO_BUILD_EXAMPLES=OFF -DVORO_BUILD_CMD_LINE=OFF -DVORO_ENABLE_DOXYGEN=OFF
cmake --build build --parallel ${nproc}
cmake --install build
install_license /workspace/srcdir/voro/LICENSE

cd ../
cd voropp-cxxwrap-julia
cmake -B build -DCMAKE_INSTALL_PREFIX=$prefix -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} -DCMAKE_INSTALL_PREFIX=$prefix -DJulia_PREFIX=$prefix -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build
"""

julia_versions = [v"1.10", v"1.11", v"1.12"]
# These are the platforms we will build for by default, unless further
# platforms are passed in on the command line
platforms = [Platform(arch, "linux"; libc = "glibc", julia_version=jv) for arch in ("i686", "x86_64", "aarch64") for jv in julia_versions]
#    Platform("i686", "linux"; libc = "glibc"),
#    Platform("x86_64", "linux"; libc = "glibc", julia_version="1.10"),
#    Platform("aarch64", "linux"; libc = "glibc"),
#    Platform("armv6l", "linux"; call_abi = "eabihf", libc = "glibc"),
#    Platform("armv7l", "linux"; call_abi = "eabihf", libc = "glibc"),
#    Platform("powerpc64le", "linux"; libc = "glibc"),
#    Platform("riscv64", "linux"; libc = "glibc"),
#    Platform("x86_64", "linux"; libc = "musl"),
#    Platform("aarch64", "linux"; libc = "musl"),
#    Platform("armv6l", "linux"; call_abi = "eabihf", libc = "musl"),
#    Platform("armv7l", "linux"; call_abi = "eabihf", libc = "musl")
#    Platform("x86_64", "windows"; )
#]

platforms = expand_cxxstring_abis(platforms)


# The products that we will ensure are always built
products = [
    LibraryProduct("libvoro++wrap", :libvoropp_wrap)
]

# Dependencies that must be installed before this package can be built
dependencies = [
    Dependency("libcxxwrap_julia_jll"; compat="0.14.3")
    BuildDependency("libjulia_jll")
]

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies
    ;
    julia_compat="1.10", preferred_gcc_version=v"9")
