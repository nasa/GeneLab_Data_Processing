FROM continuumio/anaconda3:2024.02-1

MAINTAINER Olabiyi Obayomi <obadbotanist@yahoo.com>

RUN /opt/conda/bin/conda init bash && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda update conda -y && \
    /opt/conda/bin/conda clean -afy

# copy the necessary files
COPY R.yaml ./

RUN apt-get --allow-releaseinfo-change update  && \
    apt-get upgrade -y && \
    dpkg --configure -a 

# Install environment
RUN conda env create -f R.yaml && \
    conda clean --all

RUN apt-get clean && \
    apt-get autoremove 

RUN echo "source activate /opt/conda/envs/R_env/" > ~/.bashrc

RUN apt-get install -y groff

ENV PATH="/opt/conda/envs/R_env/bin/:/opt/conda/bin:$PATH"

RUN apt-get install -y procps
CMD ["/bin/bash"]
