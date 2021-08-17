Channel
  .fromPath("$projectDir/bootstrap/meta.tsv", checkIfExists: true)
  .into { sample_meta; sample_meta_q2r }
Channel
  .fromPath("$projectDir/bootstrap/denoise/feat_tab.qza", checkIfExists: true)
  .into{ feat_tab; feat_tab_qc; feat_tab_ord }
Channel
  .fromPath("$projectDir/bootstrap/denoise/rep_seqs.qza", checkIfExists: true)
  .into{ rep_seqs; rep_seqs_pg; rep_seqs_da }
Channel
  .fromPath("$projectDir/bootstrap/phyloseq.rds", checkIfExists: true)
  .into{ ps_da; ps_abund; ps_ord; ps_ad } 
Channel
  .fromPath("http://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/V15.22/HOMD_16S_rRNA_RefSeq_V15.22.fasta")
  .set{ homd_fasta }
Channel
  .fromPath("http://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/V15.22/HOMD_16S_rRNA_RefSeq_V15.22.qiime.taxonomy")
  .set{ homd_tax }

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

process qiime_taxonomy {
    label 'qiime'

    publishDir "$baseDir/results/tax", mode: 'copy', overwrite: true

    input:
    file rep_seqs
    file homd_classifier 

    output:
    file "taxonomy.qza" into taxonomy
    val "taxonomy" into tax_channel

    """
    env TMPDIR='/data/projects/wp6/tmp' qiime feature-classifier classify-sklearn \
      --i-reads ${rep_seqs} \
      --i-classifier ${homd_classifier} \
      --p-n-jobs ${task.cpus} \
      --o-classification taxonomy.qza
    """
}

process qiime_phylogeny {
    label 'qiime'

    input:
    file rep_seqs_pg
    val tax_channel

    output:
    file "unrooted_tree.qza" into unrooted_tree
    file "rooted_tree.qza" into rooted_tree, rooted_tree_ord

    """
    env TMPDIR='/data/projects/wp6/tmp' qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences ${rep_seqs_pg} \
      --o-alignment aligned_rep_seqs.qza \
      --o-masked-alignment masked_aligned_rep_seqs.qza \
      --o-tree unrooted_tree.qza \
      --o-rooted-tree rooted_tree.qza \
      --p-n-threads ${task.cpus}
     """
}

process qiime_QC {
    label 'qiime'

    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file feat_tab_qc

    output:
    file "first_qc/"

    """
    qiime feature-table summarize \
        --i-table ${feat_tab_qc} \
        --output-dir first_qc/
    """
}

process qiime_ordination {
    publishDir "$baseDir/results", mode: 'copy', overwrite: true
    container 'qiime2/core:2019.7'

    input:
    file feat_tab_ord
    file rooted_tree_ord
    file sample_meta

    output:
    file "core-metrics-results/"

    """
    qiime diversity core-metrics-phylogenetic \
      --i-phylogeny ${rooted_tree_ord} \
      --i-table ${feat_tab_ord} \
      --p-sampling-depth 28567 \
      --m-metadata-file ${sample_meta} \
      --output-dir core-metrics-results/
    """
}

process rep_seq_to_df {
    input:
    file rep_seqs_da

    output:
    file "seq_hash.rds" into seq_hash

    """
    rep_seq_to_df.R $rep_seqs_da   
    """

}

process differential_abundance {
    container 'nebfold/ps'
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_da
    file seq_hash

    output:
    file "*.png"
    file "diff_abund.csv" into diff_abund

    """
    diff_abund.R $ps_da $seq_hash
    """
}

process plot_abundance {
    conda 'conda-forge::r-tidyverse=1.3.1 bioconda::bioconductor-phyloseq>=1.24.2 conda-forge::r-svglite>=1.2.3' 
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_abund
    file diff_abund

    output:
    file "*.png"

    """
    plt_abundance.R ${ps_abund} ${diff_abund}
    """
}

process plot_ordination {
    container 'nebfold/ps' 
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_ord

    output:
    file "*.png"
    file "*.txt"

    """
    plot_ordination.R ${ps_ord} 
    """
}

process alpha_diversity {
    conda 'conda-forge::r-tidyverse=1.3.1 bioconda::bioconductor-phyloseq>=1.24.2 conda-forge::r-svglite>=1.2.3' 
    publishDir "$baseDir/results", mode: 'copy', overwrite: true

    input:
    file ps_ad

    output:
    file "*.txt"
    file "*.png"

    """
    plot_alpha.R ${ps_ad}
    """
}
