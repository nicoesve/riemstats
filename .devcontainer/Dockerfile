FROM rocker/r-ver:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install renv first
RUN R -e "install.packages('renv')"

# Use renv to install and track other R packages
RUN R -e "renv::init(bare = TRUE); renv::install(c('here', 'testthat', 'devtools', 'roxygen2', 'pkgdown', 'tinytex')); renv::install('nicoesve/riemtan'); renv::snapshot()"

RUN R -e 'install.packages("languageserver")'

RUN R -e 'install.packages("vscDebugger", repos = "https://manuelhentschel.r-universe.dev")'

RUN apt-get update && apt-get install -y git
RUN apt-get update && apt-get install -y pandoc