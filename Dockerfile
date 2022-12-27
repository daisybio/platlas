# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.0.2
MAINTAINER Leonora Raka <leonoraraka1801@gmail.com>

# system libraries of general use
## install debian packages
RUN apt-get update && apt-get install -y \
libxml2-dev \ 
redis-server \
libcurl4-gnutls-dev \
libssl-dev \
libgmp3-dev \
libmpfr-dev \
libjpeg-dev \
libhiredis-dev \
libproj-dev \
libharfbuzz-dev \
libfribidi-dev \
libfreetype6-dev \
libpng-dev \
libtiff5-dev \
libjpeg-dev \
libc++1 \
libgfortran-5-dev \
libquadmath0 \
zlib1g \
libicu-dev \
libbz2-dev \
libresolv \
libc6-dev \
libobjc



WORKDIR /Leonora/
RUN mkdir platlas
ADD . platlas

#update shiny server conf and configure it to run hitseekr in single app mode
RUN sed -i 's/site_dir \/Leonora\/shiny-server;/app_dir \/Leonora\/Platlas;/g' /etc/shiny-server/shiny-server.conf

# go to project directory
WORKDIR /Leonora/platlas/

ENV RENV_VERSION 0.9.3-86
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org')); \
  remotes::install_github('rstudio/renv@${RENV_VERSION}'); \
  renv::restore()"
  
EXPOSE 8080

CMD ["R", "-e", "shiny::runApp('/Leonora/platlas/app.R', host = '0.0.0.0', port = 8080)"]