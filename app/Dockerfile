# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.5.0
#MAINTAINER Leonora Raka <leonoraraka1801@gmail.com>

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
libgfortran5 \
libquadmath0 \
zlib1g \
libicu-dev \
libbz2-dev \
libc6-dev \
libobjc4 \
cmake \
libmagick++-dev \
libv8-dev 


WORKDIR /Leonora/
RUN mkdir platlas
ADD . platlas

#update shiny server conf and configure it to run in single app mode
RUN sed -i 's/site_dir \/Leonora\/shiny-server;/app_dir \/Leonora\/platlas;/g' /etc/shiny-server/shiny-server.conf

# go to project directory
WORKDIR /Leonora/platlas/

ENV RENV_VERSION 1.1.5
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org')); \
  renv::restore()"
 
EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/Leonora/platlas/', host = '0.0.0.0', port = 3848)"]
