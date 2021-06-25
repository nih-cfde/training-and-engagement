This is the text that specifies software to set up in our docker container. It should be saved as "Dockerfile".

```
# Specify the base image
FROM ubuntu:20.04
MAINTAINER Marisa Lim

# Set the working directory to be used when the docker gets run
WORKDIR /usr

# disable interactive mode
ENV DEBIAN_FRONTEND noninteractive

# Do a few updates of the base system, install helper packages and R
RUN apt-get update && \
    apt-get install -y \
        unzip \
        wget \
        autoconf \
        autogen \
        make \
        g++ \
        gcc \
        git \
        automake \
        pkg-config \
        zlib1g-dev \
        curl \
        gdebi-core \
        r-base \
        r-base-dev \
        ghostscript-x

# install vcftools
RUN git clone https://github.com/vcftools/vcftools.git && \
    cd vcftools && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install

# Install the qqman library
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('qqman')"
```
