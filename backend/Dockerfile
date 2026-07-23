# One container: R serves BOTH the frontend and the analysis API.
# No separate Vercel deploy, no CORS, one URL.
#
#   docker build -t expression-console .
#   docker run -p 8000:8000 -v econsole-data:/data expression-console
#   -> open http://localhost:8000
#
# Needs >= 2 GB RAM (4 GB for large series like GSE39582 / GSE48350).
# The -v volume keeps registered accounts across restarts; without it they reset.

FROM bioconductor/bioconductor_docker:latest

# libsodium is required by the R `sodium` package (password hashing)
RUN apt-get update && apt-get install -y --no-install-recommends libsodium-dev \
 && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('plumber','jsonlite','sodium'), repos='https://cloud.r-project.org')" \
 && R -e "BiocManager::install(c('GEOquery','limma','Biobase'), update=FALSE, ask=FALSE)"

WORKDIR /app
COPY plumber.R /app/plumber.R
COPY start.R   /app/start.R
COPY www       /app/www

# accounts persist here — mount a volume to keep them
RUN mkdir -p /data
ENV USERS_FILE=/data/users.json

EXPOSE 8000
CMD ["Rscript", "/app/start.R"]
