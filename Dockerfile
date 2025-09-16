FROM rocker/r-ver:4.2.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    libpng-dev \
    libjpeg-dev \
    libmariadb-dev \
    libmariadb-dev-compat \
    && apt-get clean

# Install Bioconductor and CRAN packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(version = '3.16')"
RUN R -e "BiocManager::install(c('DESeq2', 'tidyverse', 'pheatmap'), ask=FALSE, update=TRUE)"

# Set working directory
WORKDIR /app

# Copy analysis script
COPY main.R .

# Set default command
CMD ["Rscript", "main.R"]
