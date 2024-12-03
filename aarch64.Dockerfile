# this Dockerfile is directly adapted from
# https://github.com/biocore/scikit-bio/blob/31123c6471dc62f45a55bfdff59c61a4850be367/aarch64.Dockerfile#LL1C1-L16C89
FROM --platform=linux/arm64 condaforge/linux-anvil-aarch64
RUN sudo yum update -y && \
	sudo yum install -y make git && \
	sudo yum clean all
ARG PYTHON_VERSION
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && conda create -n testing -c conda-forge --yes python=$PYTHON_VERSION gxx_linux-aarch64"
COPY . /work
WORKDIR /work
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate testing && conda env update -q -f ci/conda_host_env.yml"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate testing && conda install -q --yes --file ci/aarch64.conda_requirements.txt"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate testing && pip install -e . --no-deps"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate testing && conda list"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate testing && make test"
