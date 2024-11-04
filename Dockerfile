FROM rocker/shiny

LABEL version="1.0"
LABEL authors="t.darde@wepredic.com"
LABEL description="Docker image containing the Shiny app to visualize single cell data."


RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    git-core \
    build-essential \
    libssl-dev \
    libcurl4-gnutls-dev \
    curl \
    cmake \
    libsodium-dev \
    libxml2-dev \
    libicu-dev \
    libssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libhdf5-dev \
    default-jdk \
    libbz2-dev \
    libpcre3-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# fix for anaconda install (needed for xlsx R package)
# as per https://github.com/s-u/rJava/issues/235

RUN apt-get install -y wget 

RUN apt-get install bzip2 

RUN mkdir -p /opt/conda 
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh -O /opt/conda/miniconda.sh \ 
&& bash /opt/conda/miniconda.sh -b -p /opt/miniconda 

# Uncomment once you are ready to start productionizing the image 
# COPY environment.yaml /tmp 
# RUN . /opt/miniconda/bin/activate && conda env update --name base --file /tmp/environment.yaml 


RUN mkdir -p /usr/local/lib/R/site-library

# Définir les permissions pour le répertoire des packages R
RUN chown -R shiny:shiny /usr/local/lib/R /usr/local/lib/R/site-library

RUN R CMD javareconf
ENV R_LIBS_USER=/usr/local/lib/R/site-library
COPY ./app /shinyApp
WORKDIR /shinyApp
# RUN Rscript ./install.R

USER shiny

RUN Rscript ./installpackages.R
RUN Rscript ./reticulate_create_env.R
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

USER root
RUN chown -R shiny:shiny /shinyApp

EXPOSE 6868

CMD ["Rscript", "app.R"]