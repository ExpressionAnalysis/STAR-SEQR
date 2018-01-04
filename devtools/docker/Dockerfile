FROM ubuntu:16.04

MAINTAINER Jeff Jasper jasper1918@gmail.com

RUN apt-get -y update --fix-missing && \
    apt-get install -y curl wget gzip bzip2 unzip git \
    g++ make libboost-dev libboost-thread-dev libboost-system-dev zlib1g-dev ncurses-dev locales\
    libglib2.0-0 libxext6 libsm6 libxrender1 libxml2-dev libxslt-dev ca-certificates git gcc libdb5.3 libdb5.3-dev libcurl4-openssl-dev && \
    apt-get clean && dpkg-reconfigure locales && locale-gen en_US.UTF-8 && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

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
    rm tini.deb

# Set envinroment
ENV PATH /opt/conda/bin:$PATH
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]

# Install dependencies for STAR-SEQR
WORKDIR /opt

# Install dependencies
RUN conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda config --add channels r && \
    conda config --add channels bioconda && \
    conda install -y pip nose coverage cython pysam pandas samtools biobambam velvet ucsc-gtftogenepred salmon star gffread && \
    conda clean -ilty

# Install STAR-SEQR
ENV STARSEQR_VERSION  0.6.6
RUN wget https://github.com/ExpressionAnalysis/STAR-SEQR/archive/v${STARSEQR_VERSION}.tar.gz  && \
    tar -zxvf v${STARSEQR_VERSION}.tar.gz && \
    cd STAR-SEQR-${STARSEQR_VERSION} && \
    python setup.py build && \
    python setup.py install && \
    python setup.py clean && \
    nosetests --with-coverage --cover-package=starseqr_utils

WORKDIR /data
RUN chmod 777 /data
