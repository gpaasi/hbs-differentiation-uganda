FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('targets','yaml','purrr'), repos='https://cloud.r-project.org')"

WORKDIR /project
COPY . /project

CMD ["R", "-q", "-e", "library(targets); tar_make()"]
