# Setting up Docker

- need to make an account on Docker hub if you want to push docker there
- there's also a google storage option
- make Dockerfile (it should be called "Dockerfile")
   - basic syntax guide
- install docker
- build and test docker
- push docker to share it with others!


Installed docker and ran all commands on Mac laptop:
```
# build container
# docker build -t <user name>/<docker name>:<version number> <file location of Dockerfile>
docker build -t mlim13/gwas_test:tag0 .
# test that container works
docker run -it mlim13/gwas_test:tag0
# push to dockerhub (username: mlim13, repo: gwas_test, version: tag0)
docker push mlim13/gwas_test:tag0
```
