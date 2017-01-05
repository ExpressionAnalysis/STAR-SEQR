FROM ubuntu:14.04

MAINTAINER Jeff Jasper jasper1918@gmail.com

# Update System and packages
RUN apt-get -y update --fix-missing
RUN apt-get install -y curl wget g++ make libboost-dev libboost-thread-dev libboost-system-dev zlib1g-dev ncurses-dev unzip gzip bzip2 \
    libglib2.0-0 libxext6 libsm6 libxrender1 libxml2-dev libxslt-dev ca-certificates git gcc libdb5.1 libdb5.1-dev libcurl4-openssl-dev \
    python python-dev python-pip

# Install miniconda to /opt/conda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.1.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Use Tini
RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

ENV PATH /opt/conda/bin:$PATH

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]

# Install dependencies
WORKDIR /opt
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda install -y pip nose cython pysam==0.9.0 pandas==0.18.1 samtools biobambam star velvet ucsc-gtftogenepred

# Install STAR-SEQR
RUN wget https://github.com/ExpressionAnalysis/STAR-SEQR/archive/v0.1.0.tar.gz  && \
    tar -zxvf v0.1.0.tar.gz && \
    cd STAR-SEQR-0.1.0 && \
    python setup.py build && \
    python setup.py install && \
    python setup.py clean
