FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.1.1 as scripts

LABEL MAINTAINER="Karen Fang <karen@dnastack.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG GCLOUD_CLI_VERSION
ENV GCLOUD_CLI_VERSION "${GCLOUD_CLI_VERSION}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		build-essential \
		time \
		wget \
		gfortran \
		zlib1g-dev \
		libbz2-dev \
		liblzma-dev \
		libpcre3-dev \
		libcurl4-openssl-dev \
		libxml2-dev \
		libssl-dev \
		libfontconfig1-dev \
		libtiff-dev \
		libharfbuzz-dev \
		libfribidi-dev \
	&& rm -rf /var/lib/apt/lists/*

# R
ARG R_VERSION
ENV R_VERSION "${R_VERSION}"
RUN wget https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz && \
	tar -zxvf R-${R_VERSION}.tar.gz --directory /opt && \
	rm R-${R_VERSION}.tar.gz

ENV PATH "${PATH}:/opt/R-${R_VERSION}"

WORKDIR /opt/R-${R_VERSION}
RUN ./configure --with-readline=no --with-x=no && \
	make && \
	make install

COPY ./requirements.txt /opt/requirements.txt
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cran.r-project.org"); \
  packages <- readLines("/opt/requirements.txt"); \
  BiocManager::install(packages, ask = FALSE)'
RUN Rscript -e 'install.packages("remotes", repos="http://cran.r-project.org")'
RUN Rscript -e 'remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)'
RUN Rscript -e 'remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)'
RUN Rscript -e 'remotes::install_github("mojaveazure/seurat-disk", "seurat5", quiet = TRUE)'

COPY scripts /opt/scripts
ENV PATH "${PATH}:/opt/scripts"
