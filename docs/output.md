# Outputs

While running, Nanocompore stores the results in an SQLite database. After all transcripts are analyzed, it would perform a postprocessing step and output an easy to read TSV (tab-separated values) file.

## Results TSV

The results TSV is found in the output directory, typically named as "out_nanocompore_results.tsv". There's a set of columns that always will be present in the results and others that may appear or not, depending on the input configuration. For example, changing the statistical tests that are performed by Nanocompore would change the set of output columns.

### Mandatory output columns

| Column         | Description                                                                                                                                                                                      | Example           |
| -----------    | -------------                                                                                                                                                                                    | --------          |
| **pos**        | Starting position on the transcript's reference sequence of the k-mer for which the row reports results. Indexing is 0-based, meaning that the first nucleotide on the transcript is position 0. | 42                |
| **chr**        | Chromosome id. This will be set only if a GTF file is provided using the `gtf` parameter in the configuration.                                                                                   | chr1              |
| **genomicPos** | The position on the chromosome corresponding to the position on the transcript. This is only set when `gtf` file is provided.                                                                    | 1324576           |
| **ref_id**     | Reference name of the transcript.                                                                                                                                                                | ENST00000676788.1 |
| **strand**     | Genomic strand. This is set only when `gtf` file is provided.                                                                                                                                    | +                 |
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

## Shift statistics TSV

The shift statistics TSV gives summary statistics (mean, median, standard deviation) at the position level for the signal measurements (current intensity and dwell time) for the two conditions. The data will always be gathered during the analysis and saved to the database, but it will be exported to a TSV file only when `export_shift_stats: True` is added to the configuration. The TSV will includes the following columns:


| Column                  | Description                                                                                                                                                                                      | Example           |
| -----------             | -------------                                                                                                                                                                                    | --------          |
| **ref_id**              | Reference name of the transcript.                                                                                                                                                                | ENST00000676788.1 |
| **pos**                 | Starting position on the transcript's reference sequence of the k-mer for which the row reports results. Indexing is 0-based, meaning that the first nucleotide on the transcript is position 0. | 42                |
| **c1_mean_intensity**   | Mean value for the current intensity at the position for condition 1.                                                                                                                            | 78.91             |
| **c2_mean_intensity**   | Mean value for the current intensity at the position for condition 2.                                                                                                                            | 81.91             |
| **c1_median_intensity** | Median value for the current intensity at the position for condition 1.                                                                                                                          | 78.21             |
| **c2_median_intensity** | Median value for the current intensity at the position for condition 2.                                                                                                                          | 80.21             |
| **c1_sd_intensity**     | Standard deviation for the current intensity at the position for condition 1.                                                                                                                    | 1.21              |
| **c2_sd_intensity**     | Standard deviation for the current intensity at the position for condition 2.                                                                                                                    | 2.21              |
| **c1_mean_dwell**       | Mean value for the dwell time at the position for condition 1.                                                                                                                                   | 0.31              |
| **c2_mean_dwell**       | Mean value for the dwell time at the position for condition 2.                                                                                                                                   | 0.23              |
| **c1_median_dwell**     | Median value for the dwell time at the position for condition 1.                                                                                                                                 | 0.29              |
| **c2_median_dwell**     | Median value for the dwell time at the position for condition 2.                                                                                                                                 | 0.33              |
| **c1_sd_dwell**         | Standard deviation for the dwell time at the position for condition 1.                                                                                                                           | 0.19              |
| **c2_sd_dwell**         | Standard deviation for the dwell time at the position for condition 2.                                                                                                                           | 0.23              |

## Result database

The TSV files described above are created for the user's convenience at the end of Nanocompore's run. All data for them are sourced from the SQLite database that Nanocompore uses throughout the run to store all results. The database would be found in the output directory under the filename "out_sampComp_sql.db".

The database schema is as follows:

```sql

CREATE TABLE transcripts (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name VARCHAR NOT NULL UNIQUE
);
CREATE INDEX transcripts_name_index
    ON transcripts(name);
  
CREATE TABLE kmer_results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    transcript_id INTEGER NOT NULL,
    pos INTEGER NOT NULL,
    kmer INTEGER NOT NULL,
    c1_mean_intensity FLOAT,
    c1_mean_dwell FLOAT,
    c1_median_intensity FLOAT,
    c1_median_dwell FLOAT,
    c1_std_intensity FLOAT,
    c1_std_dwell FLOAT,
    c2_mean_intensity FLOAT,
    c2_mean_dwell FLOAT,
    c2_median_intensity FLOAT,
    c2_median_dwell FLOAT,
    c2_std_intensity FLOAT,
    c2_std_dwell FLOAT,
    UNIQUE (transcript_id, pos),
    FOREIGN KEY (transcript_id) REFERENCES transcripts(id)
);
CREATE INDEX kmer_results_transcript_id_index
    ON kmer_results(transcript_id);
```

However, depending on the choice of statistical tests, the samples used and other parameters, additional column may be added. For example, suppose we're using the GMM and KS tests for comparing 3 knock-down and 3 wilde type samples. We'd get the following additional columns in the `kmer_results` table.

```sql
CREATE TABLE kmer_results (
    ...
    -- GMM columns
    GMM_chi2_pvalue FLOAT,
    GMM_LOR VARCHAR,
    KD_1_mod FLOAT,
    KD_1_unmod FLOAT,
    KD_2_mod FLOAT,
    KD_2_unmod FLOAT,
    KD_3_mod FLOAT,
    KD_3_unmod FLOAT,
    WT_1_mod FLOAT,
    WT_1_unmod FLOAT,
    WT_2_mod FLOAT,
    WT_2_unmod FLOAT,
    WT_3_mod FLOAT,
    WT_3_unmod FLOAT,
    -- KS columns
    KS_intensity_pvalue FLOAT,
    KS_dwell_pvalue FLOAT,
    ...
)
```

