# Voxel Cores: Efficient, robust, and provably good approximation of 3D medial axes

![image](https://user-images.githubusercontent.com/3752013/204099760-80199d44-f083-4777-80b5-0738d2452060.png)

[Yajie Yan](https://yajieyan.github.io/), [David Letscher](https://cs.slu.edu/people/letscher), [Tao Ju](https://www.cse.wustl.edu/~taoju/)
<br/>*ACM Transaction on Graphics (Proceedings of SIGGRAPH 2018)*<br/>

[`Project Page`](https://yajieyan.github.io/publication/voxelcore/overview/)

The tool computes medial axis of 3d surface with union of voxels.

## Disclaimer

The repo is built & tested using `bazel` on `Mac M1` machine. Since `bazel` is cross-platform, and no platform specific symbols are used in the code, I suspect the code can be made to build on other OSes with reasonable effort, but I don't have the proper machines to test the build. I'm happy to merge your pull-requests if you have a fix.

## Repo setup

`git clone --recursive git@github.com:danielyan86129/voxel_ma.git`
This will pull the repo & all the dependencies. All deps are maintained as git submodules in `3rdparty`

## Build

Kick off the build command from within the repo folder:

`bazel build -c dbg --spawn_strategy=local //:main_voroUtility` (`-c dbg` is for building a version with debug symbols in case you want to debug the program. Otherwise feel free to omit it.)

The executable will be built under `<repo>/bazel-bin`. For details about actual arguments see [this page](https://yajieyan.github.io/publication/voxelcore/readme/)

### Why bazel
We use `bazel` to build the program, which is a modern & popular build tool created by Google. In my opinion, compared to `cmake` bazel does a better job abstracting away things that you usually don't care about, and serves most usecases well with simpler & fewer functions. See [this](https://bazel.build/install/) on how to install `bazel` for your OS. E.g. I followed https://bazel.build/install/os-x#install-on-mac-os-x-homebrew to install `bazel` on Mac.
