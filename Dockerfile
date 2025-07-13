FROM r-base:4.5.1

# Install any Linux dependencies here (optional)
RUN apt update && apt-get install -y \
libcurl4-openssl-dev \
libssl-dev \
libxml2-dev \
file \
r-cran-devtools \
bash \
bash-completion \
ncbi-blast+

ENV R_LIBS_USER=/usr/local/lib/R/site-library

RUN Rscript -e 'install.packages("pacman")'
RUN Rscript -e 'pacman::p_load(type = "binary", shiny)'
RUN Rscript -e 'pacman::p_load(type = "binary", shinyjs)'
RUN Rscript -e 'pacman::p_load(type = "binary", shinyWidgets)'
RUN Rscript -e 'pacman::p_load(type = "binary", shinyBS)'
RUN Rscript -e 'pacman::p_load(type = "binary", shinydashboard)'
RUN Rscript -e 'pacman::p_load(type = "binary", shinycssloaders)'
RUN Rscript -e 'pacman::p_load(type = "binary", DT)'
RUN Rscript -e 'pacman::p_load(type = "binary", data.table)'
RUN Rscript -e 'pacman::p_load(type = "binary", openxlsx)'
RUN Rscript -e 'pacman::p_load(type = "binary", Biostrings)'
RUN Rscript -e 'pacman::p_load(type = "binary", GenomicFeatures)'
RUN Rscript -e 'pacman::p_load(type = "binary", rtracklayer)'
RUN Rscript -e 'pacman::p_load(type = "binary", stringi)'
RUN Rscript -e 'pacman::p_load(type = "binary", stringr)'
RUN Rscript -e 'pacman::p_load(type = "binary", ggsci)'
RUN Rscript -e 'pacman::p_load(type = "binary", ggplot2)'
RUN Rscript -e 'pacman::p_load(type = "binary", ggrepel)'
RUN Rscript -e 'pacman::p_load(type = "binary", msa)'
RUN Rscript -e 'pacman::p_load(type = "binary", tidyr)'
RUN Rscript -e 'pacman::p_load(type = "binary", ggraph)'
RUN Rscript -e 'pacman::p_load(type = "binary", graphlayouts)'
RUN Rscript -e 'pacman::p_load(type = "binary", scales)'
RUN Rscript -e 'pacman::p_load(type = "binary", impute)'
RUN Rscript -e 'pacman::p_load(type = "binary", igraph)'
RUN Rscript -e 'pacman::p_load(type = "binary", scatterpie)'
RUN Rscript -e 'pacman::p_load(type = "binary", plotfunction)'
RUN Rscript -e 'pacman::p_load(type = "binary", mapplot)'
RUN Rscript -e 'pacman::p_load(type = "binary", KinSwingR)'
RUN Rscript -e 'pacman::p_load(type = "binary", gdata)'

RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install()'
RUN Rscript -e 'BiocManager::install(c("msa", "systemfonts", "Biostrings", "GenomicFeatures", "GenomicRanges", "Rsamtools", "IRanges", "rtracklayer"))'

RUN Rscript -e 'devtools::install_github("drostlab/metablastr", build_vignettes = TRUE, dependencies = TRUE)'

RUN Rscript -e 'devtools::install_github("mhahsler/rBLAST", dependencies = TRUE)'
RUN Rscript -e 'install.packages("https://github.com/wangshisheng/motifeR/raw/master/rmotifx_1.0.tar.gz", repos = NULL, type = "source")'
RUN Rscript -e 'install.packages("https://github.com/wangshisheng/motifeR/raw/master/ggseqlogo_0.1.tar.gz", repos = NULL, type = "source")' 

ENV SHELL="bash"

EXPOSE 8989

CMD ["Rscript", "-e", "devtools::load_all("/ptmore"); library(PTMoreR); PTMoreR_app()"]