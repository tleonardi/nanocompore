# Outputs

While running, Nanocompore stores the results in an SQLite database. After all transcripts are analyzed, it would perform a postprocessing step and output an easy to read TSV (tab-separated values) file.

## Results TSV

The results TSV is found in the output directory, typically named as "out_nanocompore_results.tsv". There's a set of columns that always will be present in the results and others that may appear or not, depending on the input configuration. For example, changing the statistical tests that are performed by Nanocompore would change the set of output columns.

### Mandatory output columns

| Column         | Description                                                                                                                                                                                      | Example           |
| -----------    | -------------                                                                                                                                                                                    | --------          |
| **pos**        | Starting position on the transcript's reference sequence of the k-mer for which the row reports results. Indexing is 0-based, meaning that the first nucleotide on the transcript is position 0. | 42                |
| **chr**        | Chromosome id. This will be set only if a bed file is provided using the `bed` parameter in the configuration.                                                                                   | chr1              |
| **genomicPos** | The position on the chromosome corresponding to the position on the transcript. This is only set when `bed` file is provided.                                                                    | 1324576           |
| **ref_id**     | Reference name of the transcript.                                                                                                                                                                | ENST00000676788.1 |
| **strand**     | Genomic strand. This is set only when `bed` file is provided.                                                                                                                                    | +                 |
| **ref_kmer**   | The k-mer sequence as found in the reference.                                                                                                                                                    | GCTACGT           |

### Gaussian-Mixture Models (GMM) output columns

| Column              | Description                                                                                                                                                            | Example  |
| -----------         | -------------                                                                                                                                                          | -------- |
| **GMM_chi2_pvalue** | p-value obtained from performing an association testing with a Chi-squared test between the reads' condition labels and the cluster assignments obtained from the GMM. | 0.00123  |
| **GMM_chi2_qvalue** | The p-value corrected for multiple testing using the Benjamini-Hochberg procedure.                                                                                     | 1.0      |
| **GMM_LOR** | The log odds ratio from the GMM. Suppose that we have a wildtype (WT) condition and a mod-writer knockdown (KD) condition and that we have clusters C1 and C2 detected by the GMM, then `GMM_LOR=ln((WT_C1/WT_C2)/(KD_C1/KD_C2))`, where `ln` denotes the natural logarithm. | 1.3
| **<SAMPLE\>_mod**  | The number of reads for the given k-mer assigned to the GMM cluster considered to represent the modification state.                                                    | 47       |
| **<SAMPLE\>_unmod**  | The number of reads for the given k-mer assigned to the GMM cluster considered to represent the non-modification state.                                                | 72       |


The last two columns are repeated for each of the samples provided in the input configuration YAML file. `<SAMPLE>` will be substituted with the sample label used in the configuration.

**IMPORTANT:**
It's recommended that both the q-value and the `GMM_LOR` values are used when filtering the results. The q-value provides a measurement on the probability that the separation between the two conditions is due to chance, while the LOR measures the amount of separation. As a rule of thumb, we suggest considering as modified sites with `q-value <= 0.01` **and** `|GMM_LOR| >= 0.5` (i.e. absolute value of the LOR is larger than 0.5).

## Shift stats TSV
## Result database
