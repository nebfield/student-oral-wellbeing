FROM bioconductor/release_base2:R3.6.0_Bioc3.9

RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y libcairo2-dev
RUN R -e 'install.packages(c("tidyverse", "svglite", "Hmisc"))'
RUN R -e 'BiocManager::install(c("biomformat", "phyloseq", "DESeq2"))'
