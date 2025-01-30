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
		time \
		curl \
	&& rm -rf /var/lib/apt/lists/*

# Space Ranger
ARG SPACERANGER_VERSION
ENV SPACERANGER_VERSION "${SPACERANGER_VERSION}"
# N.B. this must be passed as a build arg (using the utility script: -b SPACERANGER_DOWNLOAD_QUERY_STRING='xxyy')
#   Enter information here: https://www.10xgenomics.com/support/software/space-ranger/downloads
#   The value of SPACERANGER_DOWNLOAD_QUERY_STRING = everything following the ? in the cURL request URL (no quotes)
ARG SPACERANGER_DOWNLOAD_QUERY_STRING
RUN curl -o spaceranger-${SPACERANGER_VERSION}.tar.gz \
	"https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-${SPACERANGER_VERSION}.tar.gz?${SPACERANGER_DOWNLOAD_QUERY_STRING}" \
	&& tar -zxvf spaceranger-${SPACERANGER_VERSION}.tar.gz --directory /opt \
	&& rm spaceranger-${SPACERANGER_VERSION}.tar.gz

RUN ln -s /opt/spaceranger-${SPACERANGER_VERSION}/spaceranger /usr/local/bin/

COPY ./scripts /opt/scripts
ENV PATH "${PATH}:/opt/scripts"

COPY --from=scripts /opt/scripts /opt/scripts
