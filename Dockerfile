FROM continuumio/miniconda3:latest

LABEL maintainer="Richard Barker <richard.barker@phylo.com>"
LABEL description="Arabidopsis GWAS-Spaceflight Integration analysis environment"

WORKDIR /workspace

# Copy environment file and install dependencies
COPY environment.yml .
RUN conda env create -f environment.yml -n gwas-spaceflight && conda clean -afy

# Copy project files
COPY . /workspace/arabidopsis-gwas-spaceflight

# Activate environment
ENV PATH="/opt/conda/envs/gwas-spaceflight/bin:$PATH"

# Default command
CMD ["bash"]
