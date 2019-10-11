
first_wave_forward = file("/data/wp6/sws_microbiome/depression_resequenced_data/Ulster_redo_S0_L001_R1_001.fastq.gz")
first_wave_reverse = file("/data/wp6/sws_microbiome/depression_resequenced_data/Ulster_redo_S0_L001_R2_001.fastq.gz")
first_wave_index = file("/data/wp6/sws_microbiome/depression_resequenced_data/Ulster_redo_S0_L001_I1_001.fastq.gz")
first_wave_map = file("/data/wp6/sws_microbiome/depression_resequenced_data/mappingfile_Ulster.txt")
first_sample_meta = file("/data/wp6/sws_microbiome/first_sample_meta.tsv")

second_wave_forward = file("/data/wp6/sws_microbiome/second_wave_n30/Ulster_July2017_S0_L001_R1_001.fastq.gz")
second_wave_reverse = file("/data/wp6/sws_microbiome/second_wave_n30/Ulster_July2017_S0_L001_R2_001.fastq.gz")
second_wave_index = file("/data/wp6/sws_microbiome/second_wave_n30/Ulster_July2017_S0_L001_I1_001.fastq.gz")
second_wave_map = file("/data/wp6/sws_microbiome/second_wave_n30/mappingfile_Ulster_July2017.txt")
second_sample_meta = file("/data/wp6/sws_microbiome/second_sample_meta.tsv")

silva = file("/data/wp6/silva-132-99-nb-classifier.qza")
sws_res = file("$baseDir/bootstrap/sws.res")

process qiime_import {
    container 'qiime2/core:2019.7'

    input:
    file first_wave_forward
    file first_wave_reverse
    file first_wave_index
    file second_wave_forward
    file second_wave_reverse
    file second_wave_index

    output:
    file "first_wave.qza" into first_wave
    file "second_wave.qza" into second_wave

    """
    mkdir first_wave/
    mv $first_wave_forward first_wave/forward.fastq.gz
    mv $first_wave_reverse first_wave/reverse.fastq.gz
    mv $first_wave_index first_wave/barcodes.fastq.gz
    # use env instead of nextflow's env commands (annoying docker bug)
    env TMPDIR='/data/wp6/tmp' qiime tools import \
      --type EMPPairedEndSequences \
      --input-path first_wave/ \
      --output-path first_wave.qza
    mkdir second_wave/
    mv $second_wave_forward second_wave/forward.fastq.gz
    mv $second_wave_reverse second_wave/reverse.fastq.gz
    mv $second_wave_index second_wave/barcodes.fastq.gz
    env TMPDIR='/data/wp6/tmp' qiime tools import \
      --type EMPPairedEndSequences \
      --input-path second_wave/ \
      --output-path second_wave.qza
    """
}

process qiime_demux {
    container 'qiime2/core:2019.7'
    storeDir "/data/wp6/sws_pipeline_output/cache/"
    publishDir "/data/wp6/sws_pipeline_output/demux", mode: 'copy', overwrite: true
    beforeScript 'mkdir -p /data/wp6/tmp/'

    input:
    file first_wave
    file first_wave_map
    file second_wave
    file second_wave_map

    output:
    file "first_demuxed.qza" into demuxed_first
    file "first_demuxed_ec.qza" into demuxed_first_ec
    file "second_demuxed.qza" into demuxed_second
    file "second_demuxed_ec.qza" into demuxed_second_ec

    """
    # use env instead of nextflow's env commands (annoying docker bug)
    env TMPDIR='/data/wp6/tmp' qiime demux emp-paired \
      --i-seqs $first_wave \
      --m-barcodes-file $first_wave_map \
      --m-barcodes-column BarcodeSequence \
      --p-rev-comp-mapping-barcodes \
      --o-per-sample-sequences first_demuxed.qza \
      --o-error-correction-details first_demuxed_ec.qza
    env TMPDIR='/data/wp6/tmp' qiime demux emp-paired \
      --i-seqs $second_wave \
      --m-barcodes-file $second_wave_map \
      --m-barcodes-column BarcodeSequence \
      --p-rev-comp-mapping-barcodes \
      --o-per-sample-sequences second_demuxed.qza \
      --o-error-correction-details second_demuxed_ec.qza
    """
}

process qiime_denoise {
    container 'qiime2/core:2019.7'
    publishDir "/data/wp6/sws_pipeline_output/denoise", mode: 'copy', overwrite: true

    input:
    file demuxed_first
    file demuxed_second

    output:
    file "first_feat_tab_filt.qza" into first_feat_tab, first_feat_tab_pc, first_feat_tab_qc, first_feat_tab_ord
    file "first_rep_seqs.qza" into first_rep_seqs, first_rep_seqs_pg, first_rep_seqs_da, first_rep_seqs_pc
    file "first_stats.qza" into first_stats
    file "second_feat_tab.qza" into second_feat_tab, second_feat_tab_pc, second_feat_tab_qc
    file "second_rep_seqs.qza" into second_rep_seqs, second_rep_seqs_pg, second_rep_seqs_pc
    file "second_stats.qza" into second_stats

    """
    env TMPDIR='/data/wp6/tmp' qiime dada2 denoise-paired \
      --i-demultiplexed-seqs $demuxed_first \
      --p-trunc-len-f 240 \
      --p-trunc-len-r 225 \
      --p-n-threads 0 \
      --o-table first_feat_tab.qza \
      --o-representative-sequences first_rep_seqs.qza \
      --o-denoising-stats first_stats.qza
    env TMPDIR='/data/wp6/tmp' qiime dada2 denoise-paired \
      --i-demultiplexed-seqs $demuxed_second \
      --p-trunc-len-f 240 \
      --p-trunc-len-r 225 \
      --p-n-threads 0 \
      --o-table second_feat_tab.qza \
      --o-representative-sequences second_rep_seqs.qza \
      --o-denoising-stats second_stats.qza

    # filter from first batch (same as qza_to_ps.R)
    echo SampleID > discard_samples.tsv
    echo 238 >> discard_samples.tsv
    echo 907 >> discard_samples.tsv
    echo 763 >> discard_samples.tsv

    qiime feature-table filter-samples \
      --i-table first_feat_tab.qza \
      --m-metadata-file discard_samples.tsv \
      --p-exclude-ids \
      --o-filtered-table first_feat_tab_filt.qza
    """
}

process qiime_taxonomy {
    container 'qiime2/core:2019.7'
    publishDir "/data/wp6/sws_pipeline_output/tax", mode: 'copy', overwrite: true

    input:
    file first_rep_seqs
    file second_rep_seqs
    file silva

    output:
    file "first_taxonomy.qza" into first_taxonomy
    file "second_taxonomy.qza" into second_taxonomy
    val "taxonomy" into tax_channel

    """
    env TMPDIR='/data/wp6/tmp' qiime feature-classifier classify-sklearn \
      --i-reads $first_rep_seqs \
      --i-classifier $silva \
      --p-n-jobs 16 \
      --o-classification first_taxonomy.qza
    env TMPDIR='/data/wp6/tmp' qiime feature-classifier classify-sklearn \
      --i-reads $second_rep_seqs \
      --i-classifier $silva \
      --p-n-jobs 16 \
      --o-classification second_taxonomy.qza
    """
}

process qiime_phylogeny {
    container 'qiime2/core:2019.7'

    input:
    file first_rep_seqs_pg
    file second_rep_seqs_pg
    val tax_channel

    output:
    file "first_unrooted_tree.qza" into first_unrooted_tree
    file "first_rooted_tree.qza" into first_rooted_tree, first_rooted_tree_ord
    file "second_unrooted_tree.qza" into second_unrooted_tree
    file "second_rooted_tree.qza" into second_rooted_tree

    """
    env TMPDIR='/data/wp6/tmp' qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences $first_rep_seqs_pg \
      --o-alignment aligned_rep_seqs.qza \
      --o-masked-alignment first_masked_aligned_rep_seqs.qza \
      --o-tree first_unrooted_tree.qza \
      --o-rooted-tree first_rooted_tree.qza \
      --p-n-threads 16
    env TMPDIR='/data/wp6/tmp' qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences $second_rep_seqs_pg \
      --o-alignment second_aligned_rep_seqs.qza \
      --o-masked-alignment second_masked_aligned_rep_seqs.qza \
      --o-tree second_unrooted_tree.qza \
      --o-rooted-tree second_rooted_tree.qza \
      --p-n-threads 16
      """
}

process qiime_QC {
    publishDir "$baseDir/results", mode: 'copy', overwrite: true
    container 'qiime2/core:2019.7'

    input:
    file first_feat_tab_qc
    file second_feat_tab_qc

    output:
    file "first_qc/"
    file "second_qc/"

    """
    qiime feature-table summarize \
        --i-table $first_feat_tab_qc \
        --output-dir first_qc/
    qiime feature-table summarize \
        --i-table $second_feat_tab_qc \
        --output-dir second_qc/
    """
}

process qiime_to_phyloseq {
    container 'nebfold/bioc'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file first_feat_tab
    file first_rooted_tree
    file first_taxonomy
    file first_sample_meta

    output:
    file "ps_first.rds" into ps_first_da, ps_first_abundance, ps_first_copynum, ps_lefse

    """
    qza_to_ps.R $first_feat_tab $first_rooted_tree $first_taxonomy $first_sample_meta
    """ 
}

process qiime_to_phyloseq_val {
    container 'nebfold/bioc'

    input:
    file second_feat_tab
    file second_rooted_tree
    file second_taxonomy
    file second_sample_meta

    output:
    file "ps_second.rds" into ps_val_copynum

    """
    # TODO!
    qza_to_ps_val.R $second_feat_tab $second_rooted_tree \
      $second_taxonomy $second_sample_meta
    """
}

process qiime_ordination {
    publishDir "$baseDir/results", mode: 'copy', overwrite: true
    container 'qiime2/core:2019.7'

    input:
    file first_feat_tab_ord
    file first_rooted_tree_ord
    file first_sample_meta

    output:
    file "core-metrics-results/"

    """
    qiime diversity core-metrics-phylogenetic \
      --i-phylogeny $first_rooted_tree_ord \
      --i-table $first_feat_tab_ord \
      --p-sampling-depth 28567 \
      --m-metadata-file $first_sample_meta \
      --output-dir core-metrics-results/
    """
}

process differential_abundance {
    container 'nebfold/bioc'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_first_da
    file first_rep_seqs_da

    output:
    file "*.png"
    file "diff_abund.csv" into diff_abund

    """
    diff_abund.R $ps_first_da $first_rep_seqs_da
    """
}

process plot_abundance {
    container 'nebfold/bioc'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_first_abundance
    file diff_abund

    output:
    file "*.png"

    """
    plt_abundance.R $ps_first_abundance $diff_abund
    """
}

process adjust_copynumber {
    container 'nebfold/bioc'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_first_copynum
    file ps_val_copynum

    output:
    file "ps_copyadj.rds" into ps_copyadj, ps_copyadj_pred, ps_copyadj_ad
    file "ps_val_copyadj.rds" into ps_val_copyadj

    """
    adjust_copynumber.R $ps_first_copynum ps_copyadj.rds
    adjust_copynumber.R $ps_val_copynum ps_val_copyadj.rds
    """
}

process plot_ordination {
    container 'nebfold/bioc'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_copyadj

    output:
    file "*.png"
    file "*.txt"

    """
    plot_ordination.R $ps_copyadj
    """
}

process prediction {
    container 'nebfold/bioc'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_copyadj_pred
    file ps_val_copyadj

    """
    prediction.R $ps_copyadj_pred $ps_val_copyadj
    """
}

process qiime_picrust2 {
    container 'nebfold/picrust2'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file first_feat_tab_pc
    file first_rep_seqs_pc

    output:
    file "first_*.qza" into first_ko_picrust2
    file "biom/feature-table.biom" into mc_biom

    """
    qiime picrust2 full-pipeline \
       --i-table $first_feat_tab_pc \
       --i-seq $first_rep_seqs_pc \
       --o-ko-metagenome first_ko.qza \
       --o-ec-metagenome first_ec.qza \
       --o-pathway-abundance first_mc.qza \
       --p-threads 24 \
       --p-hsp-method mp \
       --p-max-nsti 2

    qiime tools export \
        --input-path first_mc.qza \
        --output-path biom
    """
}

process prepare_lefse {
    container 'nebfold/bioc'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_lefse
    file mc_biom

    output:
    file 'lefse.txt' into lefse_input

    """
    lefse.R $ps_lefse $mc_biom
    """
}

process lefse {
    publishDir "$baseDir/results", mode: 'copy', overwrite: true
    conda 'biobakery::lefse' 

    input:
    file lefse_input

    output:
    file 'sws.res'
    file 'sws.pdf'

    """
    format_input.py $lefse_input sws.in -c 2 -u 1 -o 1000000
    run_lefse.py sws.in sws.res -a 0.01 -w 0.01 
    plot_res.py sws.res sws.pdf --format pdf --dpi 300
    # note to self: 
    # docker run -v /path/to/workdir/:/home/linuxbrew/work -it biobakery/lefse /bin/bash
    """
}

process alpha_diversity {
    container 'nebfold/bioc'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_copyadj_ad

    output:
    file "*.png"

    """
    plot_alpha.R $ps_copyadj_ad
    """
}
