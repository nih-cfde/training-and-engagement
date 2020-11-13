# Using Docker containers in a GCP instance

**Docker containers can be loaded to a GCP VM and also to Terra. These are some basic commands. Not a step by step guide yet, still learning how to set up new containers.**

The syntax to pull in a docker image is:
```
# docker pull [image name]
# for example, this will pull in an image of the Ubuntu Linux OS
docker pull ubuntu
```

To do anything with the container:
```
# docker run [image name] [bash command]
docker run ubuntu echo "Hello World!"
```

The docker container starts up, executes command, and closes. To keep it running, start it in interactive mode:
```
docker run -it ubuntu /bin/bash
```
Now, commands are executed in the ubuntu container environment. To close, type 'exit' and you're back in the GCP VM environment. You can also keep the container open and detach from it with `Ctrl+P+Q`. To check docker container session names: `docker ps -a`.

To access the VM filesystem within the docker container, you need to mount the directory you want to use to the container:
```
# -v = specify local filesystem directory to mount
# - it = interactive mode
docker run -v ~/book:/home/book -it ubuntu /bin/bash
```

To build custom docker containers, you need a Dockerfile. Some notes on syntax:
- set the base image with 'FROM' (e.g., FROM ubuntu:20.04)
- set environment variables with 'ENV variable_name variable_value'
- run commands with 'RUN'
- connect multiple commands with '&&' and newline '\'




To build a docker container, run it, and share on dockerhub (this makes it public)
```
docker build -t <dockerhub username>/<repo name>:<container tag> <location of Dockerfile>

# for example:
# the '.' means that you're running docker commands from the folder with the Dockerfile
docker build -t mlim13/terradockertest:tag0 .
# then run container
docker run -it mlim13/terradockertest:tag0
# to share it on dockerhub (note, there's also a way to put it on Google, and make it public or private)
docker push mlim13/terradockertest:tag0
```
