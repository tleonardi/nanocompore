# Nanocompore
Software package that identifies raw signal changes between two conditions from [https://github.com/jts/nanopolish](nanopolish) resquiggled dRNA-Seq data.

## How to run
python3 nanocompore.py \
	--file1 /mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_wt/nanopolish/events.txt \
	--file2 /mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_kd/nanopolish/events.txt \
	--bam1 /mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_wt/minimap/aln_transcriptome.sorted.bam \
	--bam2 /mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_kd/minimap/aln_transcriptome.sorted.bam \
	-n 100 \
	-o out2 \
	--bedfile /mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_wt/Homo_sapiens.GRCh38.90.bed \
	--mincoverage 50
## TODO
* Add tests
* Add checks to ensure input is sorted
* Add CLI

