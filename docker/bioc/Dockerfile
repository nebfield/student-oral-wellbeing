FROM bioconductor/release_base2:R3.6.0_Bioc3.9

RUN apt-get update && apt-get upgrade -y 
RUN R -e 'install.packages(c("devtools"))'
RUN R -e 'devtools::install_github("jbisanz/qiime2R@0.99.11", dependencies = TRUE)'
