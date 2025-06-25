FROM us-central1-docker.pkg.dev/dnastack-asap-parkinsons/workflow-images/util:1.2.0 as scripts

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
		curl \
		expect \
		zlib1g-dev \
		libncursesw5-dev \
		libssl-dev \
		libsqlite3-dev \
		tk-dev \
		libgdbm-dev \
		libc6-dev \
		libbz2-dev \
		libffi-dev \
		liblzma-dev \
	&& rm -rf /var/lib/apt/lists/*

# GeoMxNGSPipeline
# Linux installation was obtained from https://university.nanostring.com/geomx-dsp-software-updates
# N.B. this must in the same directory as the Dockerfile or subdirectory relative to it
ARG GEOMX_NGS_VERSION
ENV GEOMX_NGS_VERSION "${GEOMX_NGS_VERSION}"
COPY ./GeoMxNGSPipeline_Linux_${GEOMX_NGS_VERSION}.zip /opt/GeoMxNGSPipeline_Linux_${GEOMX_NGS_VERSION}.zip
RUN cd /opt \
	&& unzip -d /opt/geomx_ngs GeoMxNGSPipeline_Linux_${GEOMX_NGS_VERSION}.zip \
	&& rm GeoMxNGSPipeline_Linux_${GEOMX_NGS_VERSION}.zip \
	&& cd /opt/geomx_ngs \
	&& tar xf GeoMxNGSPipeline.tgz
RUN chmod a+rx /opt/geomx_ngs/geomxngspipeline \
	&& touch /opt/geomx_ngs/agreement_cli_yes
ENV DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1
ENV PATH "${PATH}:/opt/geomx_ngs"
