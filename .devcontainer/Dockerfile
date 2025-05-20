FROM rocker/r-ver:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install R packages
RUN R -e "install.packages(c('Matrix', 'expm', 'purrr', 'furrr', 'R6', 'languageserver', 'here', 'testthat', 'devtools', 'roxygen2', 'pkgdown', 'tinytex'))"

RUN R -e 'install.packages("vscDebugger", repos = "https://manuelhentschel.r-universe.dev")'

RUN apt-get update && apt-get install -y git
RUN apt-get update && apt-get install -y pandoc