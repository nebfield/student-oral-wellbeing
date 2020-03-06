FROM bioconductor/release_base2:R3.6.2_Bioc3.10

RUN apt-get update && apt-get upgrade -y 
RUN R -e 'install.packages("tidyverse")'
RUN R -e 'BiocManager::install(c("biomformat", "phyloseq"))'
RUN apt-get install -y libfontconfig1-dev libcairo2-dev
RUN R -e 'install.packages("svglite")'



