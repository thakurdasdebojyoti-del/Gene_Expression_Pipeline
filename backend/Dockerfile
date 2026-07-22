# Container for the R backend (plumber API). Vercel can't run this — deploy it to
# a host that runs R: AWS App Runner / EC2, Render, Railway, or Fly.io.
#
#   docker build -t geo-de-api .
#   docker run -p 8000:8000 geo-de-api      # then point the frontend at http://localhost:8000
#
# Give it >= 2 GB RAM (4 GB for large series like GSE39582 / GSE48350).

FROM bioconductor/bioconductor_docker:latest

RUN R -e "install.packages('plumber', repos='https://cloud.r-project.org')" \
 && R -e "BiocManager::install(c('GEOquery','limma','Biobase'), update=FALSE, ask=FALSE)"

WORKDIR /app
COPY plumber.R /app/plumber.R

# Most hosts inject $PORT; default to 8000 locally.
EXPOSE 8000
CMD ["R", "-e", "plumber::pr_run(plumber::plumb('/app/plumber.R'), host='0.0.0.0', port=as.integer(Sys.getenv('PORT','8000')))"]
