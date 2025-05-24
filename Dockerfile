FROM rocker/r-ver:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install renv first
RUN R -e "install.packages('renv')"

# Use renv to install and track other R packages
RUN R -e "renv::init(bare = TRUE); renv::install(c('here', 'testthat', 'devtools', 'roxygen2', 'pkgdown', 'tinytex', 'BiocManager')); renv::install('nicoesve/riemtan')"

RUN R -e "BiocManager::install(sva)"

RUN R -e "renv::snapshot()"

RUN apt-get update && apt-get install -y git
RUN apt-get update && apt-get install -y pandoc