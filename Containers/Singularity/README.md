[TOC]

## Introduction

Singularity containers are lightweight and are based on open standards that run Linux distributions. 

Singularity was also designed around the notion of extreme mobility of computing and reproducible science. Singularity is also used to perform HPC and Cloud computing. This makes it possible to develop a research work-flow on a laboratory or a laboratory server, then bundle it to run on a departmental cluster, on a leadership class supercomputer, or in the cloud.

Singularity allows FEMPAR and its dependencies to be packaged into a standard format for identical deployment almost anywhere. FEMPAR is distributed in several ready-to-run and up-to-date Singularity containers for end-users, but we also provide the environment to compile the whole FEMPAR project for developers. Please, visit [Available images](#available-images) to see the Singularity images provided.

### Singularity installation

Singularity can only run in Linux based OS.

Note that latest Singularity major version release was ported to GO language instead of Python. Official [Singularity install docs] explain how to install Singularity in your platform.

After installing Singularity, you can check its correct funcioning by typing the following command:

```bash
$ singularity run shub://GodloveD/lolcow
```

This should download a Singularity image to current directory and show a colored cow saying a random funny message. If you can see the message, then it's correctly installed.

### Naming convention

Note that Singularity pull images to files in your local host. It does not manages a registry like Docker does. Then, naming convention only applies to Singularity registries.

A remote Singularity image name is made up of slash-separated name components, prefixed by a protocol plus a registry hostname, and finally a tag separated by a colon. Usual naming structure forms a well formed`URI`.

Currently, there are 2 official endpoints for hosting Singularity images, [Singularity-Hub](http://singularity-hub.org/) and [Singularity library](https://cloud.sylabs.io/library). Naming conventions on these platforms have some minor differences:

  - **Singularity-Hub**: `protocol://hostname/organization/repository:tag`
  - **Singularity library**: `protocol://user/project/name:[tag|@sha]`

### Pull command

Pull a remote image from a registry to a local file.

```bash
$ singularity pull [pull options...] [output file] <URI>
```

The singularity pull command will download a container from the Library (`library://`), Docker Hub (`docker://`), and also Shub (`shub://`).

This will create a single file in your host. If you don't specify the output file name,  Singularity will build it from the URI like this `organization-repository-tag`.

You can find more info in the official [Singularity pull command](https://www.sylabs.io/guides/3.2/user-guide/cloud_library.html#pulling-a-container) documentation.

### Run command

Run the user-defined default command within a container. The basic Singularity `run` command takes this form:

```bash
$ singularity run [run options...] <container>
```

You can find more info in the official [Running a container](https://www.sylabs.io/guides/3.2/user-guide/quick_start.html#running-a-container) section of the Singularity documentation.

### Exec command

Execute a custom command within a container. The basic Singularity `exec` command takes this form:

```bash
$ singularity exec [exec options...] <container> <command>
```

You can find more info in the official [Executing commands](https://www.sylabs.io/guides/3.2/user-guide/quick_start.html#executing-commands) section of the Singularity documentation.

### Hybrid MPI approach

For running *MPI parallel multinode jobs*, Singularity uses what is called *Hybrid MPI approach**. This means that an MPI distribution must be installed in both, inside the container and in the host. 

MPI must be installed inside the container because it's needed to build the parallel MPI application. MPI must be installed outside the container because it should be close to the native system to be able to spawn process natively.

```bash
$ mpirun [option] singularity exec [exec options...] <container> <command>
```

With this approach, some level of compatibility between both MPI installations is needed, in particular with MPI/PMI[x] components. Exact matching between both versions (container and host) is the most general rule to ensure that this approach works properly. 

## FEMPAR 

The quickest and port FEMPAR to an HPC is Singularity. [Singularity](https://www.sylabs.io/singularity/) is a tool designed to make it easier to create, deploy, and run applications by using containers. Singularity is HPC oriented and focused in extreme mobility and reproducible science. For more information, visit [Singularity info](https://www.sylabs.io/singularity/).

FEMPAR provides several [Singularity containers](https://cloud.sylabs.io/library/fempar) ready-to-run in your host or preferred HPC center. 


### Pull FEMPAR images

To obtain FEMPAR Singularity images, you should use [Singularity pull command](#pull-command).

```bash
$ singularity pull library://fempar/default/fempar:gnu-release_p4est-serial
```
There are several [available images](#available-images) provided by FEMPAR.

Note that FEMPAR [Release Cycle](#release-cycle) is per commit on [experimental branch](https://gitlab.com/fempar/fempar/tree/experimental/) based.

### Run FEMPAR containers

To run a FEMPAR container, you should use [Singularity run command](#run-command).

```bash
$ singularity run fempar_gnu-release_p4est-serial.sif bash
```

There are several [available images](#available-images) provided by FEMPAR.

Note that, by default, Singularity containers are not writable. If you need to persistently store data created inside the container, you need to [bind-mount](https://www.sylabs.io/guides/3.2/user-guide/bind_paths_and_mounts.html#user-defined-bind-paths) host directories into the container.

### Exec FEMPAR Tutorials

Fempar containers includes executable tests and tutorials. To execute a contained FEMPAR application, you should use [Singularity exec command](#exec-command).

```bash
$ singularity exec fempar_gnu-release_p4est-serial.sif /opt/fempar/bin/tutorial_01_steady_poisson
```

There are several [available images](#available-images) provided by FEMPAR.

Note that, by default, Singularity containers are not writable. If you need to persistently store data created inside the container, you need to [bind-mount](https://www.sylabs.io/guides/3.2/user-guide/bind_paths_and_mounts.html#user-defined-bind-paths) host directories into the container.

## Details

### Available images

FEMPAR Singularity images are hosted in *Singularity library*. [Singularity library](https://cloud.sylabs.io/library/) is a service for finding and sharing container images. Through Singularity library, FEMPAR host, and allow users to get, official pre-built and up-to-date FEMPAR Singularity images.

FEMPAR Singularity library organization is currently hosting a single project:

  - **fempar** project: contains FEMPAR library, binary tutorial andtests and all their dependencies.

FEMPAR library is delivered under several environments:

  - `fempar:gnu-debug_p4est-serial`
  - `fempar:gnu-debug_p4est-parallel`
  - `fempar:gnu-release_p4est-serial`
  - `fempar:gnu-release_p4est-parallel`

You can get more [details](#details) about these environments in [Image naming](#image-naming), [Tags](#tags) and [Components and environment](#components-and-environment) subsections.

### Release cycle

FEMPAR follows agile practices to provide new features and bugfixes as soon as possible. This concept is well known as *continuous delivery*. Continuous delivery practices try to reduce as much as possible software delivery cycles.

Continuous Delivery in FEMPAR project occurs with every new approved change in the code. This means that every commit published in [experimental branch](https://gitlab.com/fempar/fempar/tree/experimental/) of FEMPAR repository automatically creates new Singularity containers ready-to-use.

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

FEMPAR Singularity image names follow the standard naming convention (`protocol://user/project/name:tag`), see [Singularity images naming](#naming-convention), with the following components:

  - Protocol: `library`
  - User: `fempar`
  - Project: `default`
  - Name: [fempar](https://cloud.sylabs.io/library/fempar)
  - Tags: `gnu-debug_p4est-serial`, `gnu-debug_p4est-parallel`, `gnu-release_p4est-serial` and `gnu-release_p4est-parallel`

Valid FEMPAR Singularity image names are built using these components, e.g. `shub://fempar/fempar-env:gnu-release_p4est-parallel`.

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

