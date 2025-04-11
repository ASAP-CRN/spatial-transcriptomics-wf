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
		unzip \
		zip \
		time \
		wget \
		cmake \
		gfortran \
		openjdk-17-jdk \
		default-jdk \
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
		libhdf5-dev \
		hdf5-tools \
	&& rm -rf /var/lib/apt/lists/*

# R
ARG R_VERSION
ENV R_VERSION "${R_VERSION}"
RUN wget https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz \
	&& tar -zxvf R-${R_VERSION}.tar.gz --directory /opt \
	&& rm R-${R_VERSION}.tar.gz

ENV PATH "${PATH}:/opt/R-${R_VERSION}"

WORKDIR /opt/R-${R_VERSION}
RUN ./configure --with-readline=no --with-x=no \
	&& make \
	&& make install

# Must install remotes first to use install_version and install_github
RUN Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/remotes_2.5.0.tar.gz", repos=NULL, type="source")'

COPY ./requirements.txt /opt/requirements.txt
RUN Rscript -e 'lines <- readLines("/opt/requirements.txt"); \
  for (line in lines) { \
    parts <- strsplit(line, "==")[[1]]; \
    pkg <- parts[1]; \
    ver <- parts[2]; \
    remotes::install_version(pkg, version = ver, upgrade == "always", repos = "http://cran.r-project.org"); \
  }'

COPY ./bioconductor_requirements.txt /opt/bioconductor_requirements.txt
RUN Rscript -e 'lines <- readLines("/opt/bioconductor_requirements.txt"); \
  for (line in lines) { \
    parts <- strsplit(line, "==")[[1]]; \
    pkg <- parts[1]; \
    ver <- parts[2]; \
    BiocManager::install(pkg, version = ver, ask = FALSE); \
  }'

# SeuratData version 0.2.2.9001
RUN Rscript -e 'remotes::install_github("satijalab/seurat-data@3e51f44", quiet = TRUE)'
# SeuratDisk version 0.0.0.9021
RUN Rscript -e 'remotes::install_github("mojaveazure/seurat-disk@877d4e1", quiet = TRUE)'

COPY scripts /opt/scripts
ENV PATH "${PATH}:/opt/scripts"
