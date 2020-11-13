# Setting up WDL

- how to build and test a WDL locally
- basic WDL syntax


Installed Cromwell and WOMtool and ran all commands on Mac laptop:

```
wget https://github.com/broadinstitute/cromwell/releases/download/47/cromwell-47.jar
wget https://github.com/broadinstitute/cromwell/releases/download/47/womtool-47.jar
```

Also need to have Java installed: https://www.oracle.com/java/technologies/javase-jdk15-downloads.html
```
# check installed
java --version
```

test script:
```
# check WDL
java -jar womtool-47.jar validate hello.wdl

# womtool will make the basic input json file!
# but you will need to customize it
java -jar womtool-47.jar inputs hello.wdl > hello.json

# run WDL
# -i = input file
java -jar cromwell-47.jar run hello.wdl -i hello.json
```


### input json file for WDL
- this is where we define inputs
- seems that inputs have to be defined in 3 places: 1) json file, 2) workflow `call` for each task, and 3) for each task
- unclear if this file is actually necessary. upload WDL without json to Terra.


### outputs
- each task gets a call-<task name> directory that has `execution` and `inputs` directories
   - the outputs, logs, stderr files are in `execution`
   - when defining the inputs in WDL script, check how the inputs are partitioned in the `inputs` directory. this was the issue for plink commands and the reason i had to define each input file flag for plink.
 
