How to get and use NEOS?
====

Developper Documentation
------------------------

[Documentation is Here](https://memphis.gitlabpages.inria.fr/neos/)

Get Neos
--------
To use last development state of Neos, please clone the master branch.

To get sources please use these commands:

      git clone git@gitlab.inria.fr:memphis/neos.git
      cd neos
      git submodule init
      git submodule update


Use Neos
-----------------------

Neos is available as a Docker Container.

- On your Operating System, install [Docker Desktop](https://docs.docker.com/get-docker/).

- Please refer to the documentation in the README at https://gitlab.inria.fr/memphis/neos for instructions on how to run Docker.

If Docker is not functioning properly or is utilizing excessive volume capacity, consider using Podman as an alternative:

   * In Linux : sudo apt-get install podman
   * In macOS : brew install podman

- To test your work on Podman, such as a file named 'code.cpp' or a folder titled 'code' located on your local machine:

```
$ podman --storage-opt ignore_chown_errors=true run -v /path/towards/your/file/folder:/folder -u root -it registry.gitlab.inria.fr/memphis/neos bash
```

- After using the command mentioned above, you will be placed in a virtual machine where you can execute and run your work.

- The next steps are to use the command 'cmake .' followed by 'make' to execute your code.

