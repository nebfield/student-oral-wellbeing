FROM bioconductor/release_base2:R3.6.0_Bioc3.9

RUN apt-get update && apt-get upgrade -y 
RUN R -e 'install.packages(c("devtools", "tidyverse"))'
RUN R -e 'BiocManager::install(c("biomformat", "yaml", "phyloseq", "Hmisc"))'
RUN R -e 'devtools::install_github("jbisanz/qiime2R")'

# TODO: tidy up below
RUN R -e 'BiocManager::install(c("DESeq2", "picante"))'
RUN R -e 'devtools::install_github("r-lib/lifecycle")' # manual install fixes weird bug
RUN apt-get install -y libpq-dev
RUN R -e 'devtools::install_github(c("tidyverse/dplyr", "tidyverse/tidyr"), dependencies = T)'
RUN R -e 'BiocManager::install("ALDEx2")'