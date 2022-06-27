# Install BLAST

## Step 1: Launch AWS instance

If you do not have access to a command line terminal on your computer (e.g., WindowsOS), you can set up an AWS remote instance to run command line programs.

Boot a `t2.micro` instance on AWS and connect your shell prompt.

!!! Tip

	If you need assistance setting up an instance check out the tutorial [Intro to AWS](../../Cloud-Platforms/Introduction-to-AWS/index.md)!


## Step 2: Install BLAST software using command line
Now, install the NCBI BLAST software. Copy and paste the following installation commands:

```
sudo apt-get update && sudo apt-get -y install python ncbi-blast+
```

The `sudo apt-get update` command is used to download package information from all configured sources. So when you run the `update` command, it downloads the package information from the internet. This updates the software list and installs the Python programming language and NCBI BLAST+.

Check that the installation was successful:

=== "Input"

	```
	blastp -version
	```

=== "Expected Output"

	```
	ubuntu@ip-172-31-17-217:~$ blastp -version
	blastp: 2.9.0+
 	Package: blast 2.9.0, build Sep 30 2019 01:57:31
	```
