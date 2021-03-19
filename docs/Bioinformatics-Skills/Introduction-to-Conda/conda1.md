# Why Conda?

![](https://i.imgur.com/2R6UDu3.jpg)

- OS-agnostic tool for package and environment management (meaning it works for Windows, Mac, and Linux!)
- clean software installation along with dependencies
- searches for compatible software versions
- many packages available - python, R, etc.
- user has admin permissions, good for use on HPC where you'd otherwise need to request software installations


### What is an Environment?

Think of software environments like school classrooms, where the school is your laptop OS. Schools have many classrooms and while they are independent rooms they are still connected by hallways. Students can walk into any room and each room will have its own set of supplies (i.e., desks, chairs, whiteboard, decorations, etc.).

For conda, the equivalent of a classroom and its supplies is a conda environment on your computer and its set of specifically downloaded software/tools. Input/output files, like students, can easily move between environments and just like a school can have multiple classrooms, conda can manage multiple conda environments. This is what makes it a great tool for maintaining different versions of the same software all on 1 computer!

### Why do we need isolated software environments?

- **versioning**: your research/project requires specific versions of software packages
- **reproducibility**: you want to experiment with new packages or features of existing packages without compromising your current workflow
- **repeatability**: you want to replicate results from a publication
- **compatibility**: your collaborators require you to use certain software/tools that are incompatibile with your current system
