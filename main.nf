
first_wave_forward = file("/data/projects/wp6/sws_microbiome/depression_resequenced_data/Ulster_redo_S0_L001_R1_001.fastq.gz")
first_wave_reverse = file("/data/projects/wp6/sws_microbiome/depression_resequenced_data/Ulster_redo_S0_L001_R2_001.fastq.gz")
first_wave_index = file("/data/projects/wp6/sws_microbiome/depression_resequenced_data/Ulster_redo_S0_L001_I1_001.fastq.gz")
first_wave_map = file("/data/projects/wp6/sws_microbiome/depression_resequenced_data/mappingfile_Ulster.txt")
first_sample_meta = file("/data/projects/wp6/sws_microbiome/first_sample_meta.tsv")

silva = file("/data/projects/wp6/silva-132-99-nb-classifier.qza")
homd_fasta = Channel.fromPath("http://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/V15.22/HOMD_16S_rRNA_RefSeq_V15.22.fasta")
homd_tax = Channel.fromPath("http://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/V15.22/HOMD_16S_rRNA_RefSeq_V15.22.qiime.taxonomy")

process qiime_homd {
    label 'qiime'
    stageInMode 'copy'

    input:
    file homd_fasta
    file homd_tax

    output:
    file "ref-seqs.qza" into homd_ref_seqs
    file "ref-taxonomy.qza" into homd_ref_tax

    """
    # add Z to all IDs
    # yes, seriously 
    sed 's/^>/>Z/' HOMD_16S_rRNA_RefSeq_V15.22.fasta > homd_z.fasta
    sed 's/^/Z/' HOMD_16S_rRNA_RefSeq_V15.22.qiime.taxonomy > tax.txt
    qiime tools import \
     --type 'FeatureData[Sequence]' \
     --input-path homd_z.fasta \
     --output-path homd.qza

    # 314F / 806 R primers
    # don't set --p-trunc-len with paired-end sequencing
    qiime feature-classifier extract-reads \
      --i-sequences homd.qza \
      --p-f-primer CCTACGGGAGGCAGCAG \
      --p-r-primer GGACTACHVGGGTWTCTAAT \
      --o-reads ref-seqs.qza

    qiime tools import \
      --type 'FeatureData[Taxonomy]' \
      --input-format HeaderlessTSVTaxonomyFormat \
      --input-path tax.txt \
      --output-path ref-taxonomy.qza
    """
}

process qiime_train {
    label 'qiime'

    input:
    file homd_ref_seqs
    file homd_ref_tax

    output:
    file "classifier.qza" into homd_classifier

    """
    qiime feature-classifier fit-classifier-naive-bayes \
      --i-reference-reads $homd_ref_seqs \
      --i-reference-taxonomy $homd_ref_tax \
      --o-classifier classifier.qza 
    """

}

process qiime_import {
    label 'qiime'

    input:
    file first_wave_forward
    file first_wave_reverse
    file first_wave_index

    output:
    file "first_wave.qza" into first_wave

    """
    mkdir first_wave/
    mv $first_wave_forward first_wave/forward.fastq.gz
    mv $first_wave_reverse first_wave/reverse.fastq.gz
    mv $first_wave_index first_wave/barcodes.fastq.gz
    # use env instead of nextflow's env commands (annoying docker bug)
    env TMPDIR='/data/projects/wp6/tmp' qiime tools import \
      --type EMPPairedEndSequences \
      --input-path first_wave/ \
      --output-path first_wave.qza
   """
}

process qiime_demux {
    label 'qiime'

    storeDir "/data/projects/wp6/sws_pipeline_v2/cache/"
    publishDir "$baseDir/results/demux", mode: 'copy', overwrite: true
    beforeScript 'mkdir -p /data/projects/wp6/tmp/'

    input:
    file first_wave
    file first_wave_map

    output:
    file "first_demuxed.qza" into demuxed_first
    file "first_demuxed_ec.qza" into demuxed_first_ec

    """
    # use env instead of nextflow's env commands (annoying docker bug)
    env TMPDIR='/data/projects/wp6/tmp' qiime demux emp-paired \
      --i-seqs $first_wave \
      --m-barcodes-file $first_wave_map \
      --m-barcodes-column BarcodeSequence \
      --p-rev-comp-mapping-barcodes \
      --o-per-sample-sequences first_demuxed.qza \
      --o-error-correction-details first_demuxed_ec.qza
    """
}

process qiime_denoise {
    label 'qiime'

    publishDir "$baseDir/results/denoise", mode: 'copy', overwrite: true

    input:
    file demuxed_first

    output:
    file "first_feat_tab_filt.qza" into first_feat_tab, first_feat_tab_pc, first_feat_tab_qc, first_feat_tab_ord, first_feat_tab_pred
    file "first_rep_seqs.qza" into first_rep_seqs, first_rep_seqs_pg, first_rep_seqs_da, first_rep_seqs_pc
    file "first_stats.qza" into first_stats

    """
    env TMPDIR='/data/projects/wp6/tmp' qiime dada2 denoise-paired \
      --i-demultiplexed-seqs $demuxed_first \
      --p-trunc-len-f 240 \
      --p-trunc-len-r 225 \
      --p-n-threads 0 \
      --o-table first_feat_tab.qza \
      --o-representative-sequences first_rep_seqs.qza \
      --o-denoising-stats first_stats.qza

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
    label 'qiime'

    publishDir "$baseDir/results/tax", mode: 'copy', overwrite: true

    input:
    file first_rep_seqs
    file homd_classifier 

    output:
    file "first_taxonomy.qza" into first_taxonomy
    val "taxonomy" into tax_channel

    """
    env TMPDIR='/data/projects/wp6/tmp' qiime feature-classifier classify-sklearn \
      --i-reads $first_rep_seqs \
      --i-classifier ${homd_classifier} \
      --p-n-jobs 16 \
      --o-classification first_taxonomy.qza
    """
}

process qiime_phylogeny {
    label 'qiime'

    input:
    file first_rep_seqs_pg
    val tax_channel

    output:
    file "first_unrooted_tree.qza" into first_unrooted_tree
    file "first_rooted_tree.qza" into first_rooted_tree, first_rooted_tree_ord

    """
    env TMPDIR='/data/projects/wp6/tmp' qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences $first_rep_seqs_pg \
      --o-alignment aligned_rep_seqs.qza \
      --o-masked-alignment first_masked_aligned_rep_seqs.qza \
      --o-tree first_unrooted_tree.qza \
      --o-rooted-tree first_rooted_tree.qza \
      --p-n-threads 16
     """
}

process qiime_QC {
    label 'qiime'

    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file first_feat_tab_qc

    output:
    file "first_qc/"

    """
    qiime feature-table summarize \
        --i-table $first_feat_tab_qc \
        --output-dir first_qc/
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
    file "ps_first.rds" into ps_first_da, ps_first_abundance, ps_lefse, ps_first_ord, ps_first_pred, ps_first_ad
    file "samp_dat.tsv" into first_samp_dat

    """
    qza_to_ps.R $first_feat_tab $first_rooted_tree $first_taxonomy $first_sample_meta
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

process rep_seq_to_df {
    container 'nebfold/bioc'

    input:
    file first_rep_seqs_da

    output:
    file "seq_hash.rds" into seq_hash

    """
    rep_seq_to_df.R $first_rep_seqs_da   
    """

}

process differential_abundance {
    container 'nebfold/ps'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_first_da
    file seq_hash

    output:
    file "*.png"
    file "*.svg"
    file "diff_abund.csv" into diff_abund

    """
    diff_abund.R $ps_first_da $seq_hash
    """
}

process plot_abundance {
    container 'nebfold/tidyr' 
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_first_abundance
    file diff_abund

    output:
    file "*.png"
    file "*.svg"

    """
    plt_abundance.R $ps_first_abundance $diff_abund
    """
}

process plot_ordination {
    container 'nebfold/ps' 
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_first_ord

    output:
    file "*.png"
    file "*.svg"
    file "*.txt"

    """
    plot_ordination.R $ps_first_ord
    """
}

process alpha_diversity {
    conda 'r::r-tidyverse=1.2.1 bioconda::bioconductor-phyloseq=1.24.2 conda-forge::r-svglite=1.2.3' 
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_first_ad

    output:
    file "*.svg"
    file "*.txt"
    file "*.png"

    """
    plot_alpha.R $ps_first_ad
    """
}
