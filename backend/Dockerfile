# Container for the R backend (plumber API). Runs on any host that runs Docker
# (Render / Railway / Fly / AWS). Vercel can't run this.
#
# Give it >= 2 GB RAM for large series (GSE39582 / GSE48350); 512 MB handles the
# small ones (GSE10072, GSE55235).

FROM bioconductor/bioconductor_docker:latest

RUN R -e "install.packages('plumber', repos='https://cloud.r-project.org')" \
 && R -e "BiocManager::install(c('GEOquery','limma','Biobase'), update=FALSE, ask=FALSE)"

WORKDIR /app
COPY plumber.R /app/plumber.R
COPY start.R   /app/start.R

# Render injects $PORT; start.R binds to it. 8000 is only the local default.
EXPOSE 8000
CMD ["Rscript", "/app/start.R"]
