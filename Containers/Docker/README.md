## Introduction

Docker containers are lightweight and are based on open standards that run on all major Linux distributions, macOS and Microsoft Windows platforms. For more information, visit [wath-docker](https://www.docker.com/what-docker).

Docker allows FEMPAR and its dependencies to be packaged into a standard format for identical deployment almost anywhere. FEMPAR is distributed in several ready-to-run and up-to-date Docker containers for end-users, but we also provide the environment to compile the whole FEMPAR project for developers. Please, visit [Available images](#available-images) to see the Docker images provided.

### Docker installation

Docker can run in Linux, Mac and Windows. Please take a look to [supported-platforms](https://docs.docker.com/install/#supported-platforms) to check if your system is supported.

Instructions to install Docker in your platform can be read in the following links:

 - Linux: [CentOS](https://docs.docker.com/install/linux/docker-ce/centos/), [Debian](https://docs.docker.com/install/linux/docker-ce/debian/), [Fedora](https://docs.docker.com/install/linux/docker-ce/fedora/), [Ubuntu](https://docs.docker.com/install/linux/docker-ce/ubuntu/)
 - Mac: [MacOS](https://docs.docker.com/docker-for-mac/install/)
 - Windows: [Windows 10](https://docs.docker.com/docker-for-windows/install/)
 
After installing Docker, you can check its correct funcioning by typing the following command:

```bash
$ docker run hello-world
```

You should see a message saying that your Docker installation is working correctly. You can get more details on how to [test Docker installation](https://docs.docker.com/get-started/#test-docker-installation).

### Naming convention

A Docker image name is made up of slash-separated name components, optionally prefixed by a registry hostname, and finally a tag separated by a colon. Usual naming structure contains the following components: `hostname/organization/repository:tag`

### Pull command

Pull an image or a repository from a registry.

```bash
$ docker pull [OPTIONS] NAME[:TAG|@DIGEST]
```

You can find more info in the official [Docker pull command](https://docs.docker.com/engine/reference/commandline/pull/) documentation.

### Run command

The basic Docker `run` command takes this form:

```bash
$ docker run [OPTIONS] IMAGE[:TAG|@DIGEST] [COMMAND] [ARG...]
```

The Docker run command must specify an IMAGE to derive the container from. 

You can find more info in the official [Docker run command](https://docs.docker.com/engine/reference/run/) documentation.

## FEMPAR 

The quickest and easiest way to start with FEMPAR is using [Docker](https://www.docker.com). [Docker](https://opensource.com/resources/what-docker) is a tool designed to make it easier to create, deploy, and run applications by using containers.

FEMPAR provides a Docker [container with the required environment](https://hub.docker.com/u/fempar) to compile the project source code and to run tutorials and tests. 


### Pull FEMPAR images

To obtain FEMPAR Docker images, you should use [Docker pull command](#pull-command).

```bash
$ docker pull fempar/fempar-env:gnu-debug_p4est-serial
```

There are several [available images](#available-images) provided by FEMPAR.

To check if your local image is up-to-date with the one hosted on Docker Hub, you can run the `pull` command. If your image is up-to-date you will get the following message, if not the image will be pulled again:

```bash
gnu-debug_p4est-serial: Pulling from fempar/fempar-env
Digest: sha256:xxxxxxxxxxxxxxxxxxxxxxxx
Status: Image is up to date for fempar/fempar-env:gnu-debug_p4est-serial
```

Note that FEMPAR [Release Cycle](#release-cycle) is per commit on [experimental branch](https://gitlab.com/fempar/fempar/tree/experimental/) based.

### Run FEMPAR containers

To run FEMPAR Docker container, you should use [Docker run command](#run-command).

```bash
$ docker run -ti fempar/fempar-env:gnu_debug_p4est-serial
```

There are several [available images](#available-images) provided by FEMPAR.

### Compile FEMPAR

Please, follow the steps below to compile FEMPAR using the Docker container:

```bash
$ WORKDIR=/data
$ SOURCES_DIR=$WORKDIR/sources
$ FEMPAR_DIR=$WORKDIR/FEMPAR
$ git clone --recursive https://github.com/fempar/fempar $SOURCES_DIR
$ cd $WORKDIR
# invokes CMake while setting up appropriate values for CMake variables
# run $SOURCES_DIR/Tools/configure -h to get an informative message on screen
$ $SOURCES_DIR/Tools/configure -s $SOURCES_DIR/SuperBuild -c GNU --without-tests 
$ make
```

### Compile FEMPAR Tutorials

If you are new to FEMPAR, the very first point to start with are the [tutorial programs](https://gitlab.com/fempar/fempar/tree/experimental/Tutorials) available at the official FEMPAR repository. 
After completing the previous compilation steps, you can compile and run FEMPAR tutorials with the following steps:

```bash
...
$ FEMPAR_TUTORIALS_DIR=$SOURCES_DIR/Tutorials
$ mkdir -p $FEMPAR_TUTORIALS_DIR
$ cd $FEMPAR_TUTORIALS_DIR
$ cmake -DFEMPAR_DIR=$FEMPAR_DIR -DFEMPAR_TUTORIAL=tutorial_01_poisson_sharp_circular_wave $SOURCES_DIR/Tutorials
$ make -j 4
$ bin/tutorial_01_poisson_sharp_circular_wave --help                              # get informative message on screen
$ bin/tutorial_01_poisson_sharp_circular_wave [optional command line arguments]   # execute the tutorial
```

This particular set of commands compiles the tutorial program named `tutorial_01_poisson_sharp_circular_wave`. You may use any of the tutorial names at `$SOURCES_DIR/Tutorials` as well.

At present, we only offer a reduced set of tutorial programs which show the usage of the most simple FEMPAR features. However, we plan in the near future to extend the current 
tutorial suite towards demonstration of the various aspects of the library. 

In the meantime, you can also take a look at the 
[serial](https://gitlab.com/fempar/fempar/tree/experimental/Sources/Tests/Serial) and [MPI-parallel test programs](https://gitlab.com/fempar/fempar/tree/experimental/Sources/Tests/Par) available at the official FEMPAR repository.
While these programs go far beyond the current tutorial programs in exploiting many of the most advanced FEMPAR features, they are, though, not fully documented,
so they are only recommended for advanced users.  Test programs are compiled as a final stage of the compilation of FEMPAR. In order to activate the compilation of tests, you have to replace
`--without-tests` by `--with-tests` in the steps above.

## Details

### Available images

FEMPAR Docker images are hosted in Docker Hub. [Docker Hub](https://hub.docker.com/) is a service provided by Docker for finding and sharing container images. Through Docker Hub, FEMPAR host, and allow users to get, official pre-built and up-to-date FEMPAR Docker images.

FEMPAR Docker Hub organization is currently hosting two available repositories:

  - **fempar** repository: contains FEMPAR library and all its dependencies
  - **fempar-env** repository: contains all required dependencies to compile FEMPAR.

#### Library

[![fempar](http://docker-badges.webbedlam.com/image/fempar/fempar)](https://hub.docker.com/r/fempar/fempar)

FEMPAR library is delivered under several environments:

  - `fempar/fempar:gnu-debug_p4est-serial`
  - `fempar/fempar:gnu-debug_p4est-parallel`
  - `fempar/fempar:gnu-release_p4est-serial`
  - `fempar/fempar:gnu-release_p4est-parallel`

You can get more [details](#details) about these environments in [Image naming](#image-naming), [Tags](#tags) and [Components and environment](#components-and-environment) subsections.

#### Development environment

[![fempar-env](http://docker-badges.webbedlam.com/image/fempar/fempar-env)](https://hub.docker.com/r/fempar/fempar-env)

FEMPAR compilation environment is delivered under several environments:

  - `fempar/fempar-env:gnu-debug_p4est-serial`
  - `fempar/fempar-env:gnu-debug_p4est-parallel`
  - `fempar/fempar-env:gnu-release_p4est-serial`
  - `fempar/fempar-env:gnu-release_p4est-parallel`    

You can get more [details](#details) about these environments in [Image naming](#image-naming), [Tags](#tags) and [Components and environment](#components-and-environment) subsections.

### Release cycle

FEMPAR follows agile practices to provide new features and bugfixes as soon as possible. This concept is well known as *continuous delivery*. Continuous delivery practices try to reduce as much as possible software delivery cycles.

Continuous Delivery in FEMPAR project occurs with every new approved change in the code. This means that every commit published in [experimental branch](https://gitlab.com/fempar/fempar/tree/experimental/) of FEMPAR repository automatically creates new Docker containers ready-to-use.

### Components and environment

All images are built on top of the official [GCC containers](https://hub.docker.com/_/gcc/). In addition all of them contain the following software components:

| Component                                                     |     Version      |
|:-------------------------------------------------------------:|:----------------:|
| [GCC](https://gcc.gnu.org/)                                   | 9.1              |
| [MKL](https://software.intel.com/en-us/mkl)                   | 64bit-2019.3-062 |
| [OpenMPI](https://www.open-mpi.org/)                          | 3.1.3            |
| [json-fortran](https://github.com/jacobwilliams/json-fortran) | 7.0.0            |
| [PETSC](https://www.mcs.anl.gov/petsc/)                       | 3.9.2            |
| [QHULL](http://www.qhull.org/)                                | 2015-src-7.2.0   |
| [HDF5](https://support.hdfgroup.org/HDF5/)                    | 1.8.21           |
| [P4EST](http://www.p4est.org/)                                | 2.2              |

### Image naming

FEMPAR Docker image names follow the standard naming convention (`hostname/organization/repository:tag`), see [Docker images naming](#naming-convention), with the following components:

  - Hostname: `registry-1.docker.io`, can be skipped a it's thee default one.
  - Organization: [fempar](https://hub.docker.com/r/fempar)
  - Repositories: [fempar](https://hub.docker.com/r/fempar/fempar) and [fempar-env](https://hub.docker.com/r/fempar/fempar-env)
  - Tags: `gnu-debug_p4est-serial`, `gnu-debug_p4est-parallel`, `gnu-release_p4est-serial` and `gnu-release_p4est-parallel`

Valid FEMPAR Docker image names are built using these components, e.g. `fempar/fempar-env:gnu-release_p4est-parallel`.

Subsection [Tags](#tags) provides more lights on tags meaning.

#### Tags

FEMPAR tags describe the environment and software installed within the Docker images. Currently there are 4 tags available:

  - `gnu-debug_p4est-serial`
  - `gnu-debug_p4est-parallel`
  - `gnu-release_p4est-serial`
  - `gnu-release_p4est-parallel`    

There are several key words to understand the meaning of this tags:

 - **debug**: Tags containing `debug` specify that all the included FEMPAR dependencies are compiled in `DEBUG` mode. 
 - **release**: Tags containing `release` specify that all the included FEMPAR dependencies are compiled in `RELEASE` mode. 
 - **p4est_serial**: Tags containing `p4est_serial` specify that `P4EST` library was compiled without MPI support. 
 - **p4est_parallel**: Tags containing `p4est_parallel` specify that `P4EST` library was compiled withMPI support. 

