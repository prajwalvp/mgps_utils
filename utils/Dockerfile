FROM ubuntu:20.04

# Install prerequisites
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq && \
    apt-get -y --no-install-recommends install \
    autotools-dev \
    autoconf \
    automake \
    build-essential \
    g++ \
    gfortran \
    git \
    python3-dev \
    python3-pip \
    rsync \
    ssh \
    vim \
    wget && \
    apt-get clean all && \
    rm -r /var/lib/apt/lists/*


# Install python dependancies
RUN pip install numpy \
    scipy \
    matplotlib \
    astropy \
    astroquery \
    requests \
    path \
    pandas \
    psrqpy \
    requests_html \
    paramiko \
    jupyter  

ENV PSRHOME /software/

#PSRCAT
WORKDIR $PSRHOME
RUN wget http://www.atnf.csiro.au/people/pulsar/psrcat/downloads/psrcat_pkg.tar.gz && \
    tar -xvf psrcat_pkg.tar.gz -C $PSRHOME
ENV PSRCAT_FILE $PSRHOME/psrcat_tar/psrcat.db
ENV PATH $PATH:$PSRHOME/psrcat_tar
WORKDIR $PSRHOME/psrcat_tar
RUN /bin/bash makeit && \
    rm -f ../psrcat_pkg.tar.gz


    
# Clone mgps repository
WORKDIR $PSRHOME
RUN git clone https://github.com/prajwalvp/mgps_utils.git
WORKDIR $PSRHOME/mgps_utils
RUN git checkout apsuse

# Adjust env paths
ENV MGPS_UTILS $PSRHOME/mgps_utils
ENV PYTHONPATH $PYTHONPATH:$MGPS_UTILS/utils
ENV PATH $PATH:$MGPS_UTILS/utils

# Make scripts executable
WORKDIR $MGPS_UTILS/utils
RUN chmod +x get_psrs_in_field.py metaxtractor.py candy_scraper2.py unpack_merge.py    
