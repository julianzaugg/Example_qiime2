# Example script to run QIIME2 on reads and generate basic results for post analysis. 
# Currently assumes forward and reverse reads.
# Assumes CasavaOneEightSingleLanePerSampleDirFmt for read filenames, i.e. SampleName_S1_L001_R1_001.fastq.gz, SampleName_S1_L001_R2_001.fastq.gz, etc...

module load miniconda3/1.1
conda activate qiime2-2020.11

READ_DIR="my/read/dir" # Location of reads, change as required.
QIIME_DIR="my/qiime/dir" # Output directory for generated results
mkdir -p $QIIME_DIR

SILVA_SEQ_DB="/srv/sw/acepipe_config/qiime_db/silva_seq_db_r138.qza"
SILVA_TAXONOMY_DB="/srv/sw/acepipe_config/qiime_db/silva_taxonomy_db_r138.qza"

ID=my_project # Prefix ID, change as required.
THREADS=36 # Number of threads/CPUs

FORWARD_TRUNC=270
REVERSE_TRUNC=260

# ------------------------------------------------------------------------
# Import reads into QIIME. Can either import directory or by constructing a manifest file

# Build manifest file. Assume that sample files are "_" delimited with the sample ID as the first part of the filename.
if [[ ! -f $QIIME_DIR/${ID}_sequences.qza ]]; then
echo "sample-id,absolute-filepath,direction" > $QIIME_DIR/manifest.csv
ls -1 $READ_DIR/*_R1_*gz | sed 's!.*/!!' | awk -v read_dir=$READ_DIR '{split($0, sample, "_"); print sample[1]","read_dir"/"$0",forward"}' >> $QIIME_DIR/manifest.csv
ls -1 $READ_DIR/*_R2_*gz | sed 's!.*/!!' | awk -v read_dir=$READ_DIR '{split($0, sample, "_"); print sample[1]","read_dir"/"$0",reverse"}' >> $QIIME_DIR/manifest.csv

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $QIIME_DIR/manifest.csv \
  --output-path $QIIME_DIR/${ID}_sequences.qza \
  --input-format PairedEndFastqManifestPhred33

# Alternatively, import reads directly. Assumes Casava format (standard Illumina namining scheme)
#qiime tools import --type SampleData[PairedEndSequencesWithQuality] \
#--input-path $READ_DIR \
#--input-format CasavaOneEightSingleLanePerSampleDirFmt \
#--output-path $QIIME_DIR/${ID}_sequences.qza
fi

# ------------------------------------------------------------------------
# Create a sample metadata file for QIIME. Simply the sample ID
if [[ ! -f $QIIME_DIR/${ID}_sample_metadata.tsv ]]; then
echo "Sample ID" > $QIIME_DIR/${ID}_sample_metadata.tsv
ls -1 $READ_DIR/*_R1_* | sed 's!.*/!!' | awk '{split($0, sample, "_"); print sample[1]}' >> $QIIME_DIR/${ID}_sample_metadata.tsv
fi
# ------------------------------------------
# Denoise and determine features

#--p-trim-left-f 10 \
#--p-trim-left-r 10 \

if [[ ! -f $QIIME_DIR/${ID}_feature_table.qza ]]; then

qiime dada2 denoise-paired \
--i-demultiplexed-seqs $QIIME_DIR/${ID}_sequences.qza \
--o-table $QIIME_DIR/${ID}_feature_table.qza \
--o-representative-sequences $QIIME_DIR/${ID}_feature_sequences.qza \
--p-trunc-len-f $FORWARD_TRUNC \
--p-trunc-len-r $REVERSE_TRUNC \
--p-n-threads $THREADS \
--o-denoising-stats $QIIME_DIR/${ID}_dada2_stats.qza

fi

# Generate DADA2 summary stats
if [[ ! -f $QIIME_DIR/${ID}_dada2_stats.tsv ]]; then

#mkdir -p $QIIME_DIR/dada2_stats
qiime tools export \
--input-path $QIIME_DIR/${ID}_dada2_stats.qza \
--output-path $QIIME_DIR/dada2_stats
mv $QIIME_DIR/dada2_stats/stats.tsv $QIIME_DIR/${ID}_dada2_stats.tsv
rm -r $QIIME_DIR/dada2_stats

fi

# Export the feature sequence summary table
if [[ ! -d $QIIME_DIR/${ID}_qiime_feature_table ]]; then

qiime tools export \
--input-path $QIIME_DIR/${ID}_feature_table.qza \
--output-path $QIIME_DIR/${ID}_qiime_feature_table

fi

# Generate the relative feature-table
if [[ ! -f $QIIME_DIR/${ID}_relative_feature_table.qza ]]; then

qiime feature-table relative-frequency \
--i-table $QIIME_DIR/${ID}_feature_table.qza \
--o-relative-frequency-table $QIIME_DIR/${ID}_relative_feature_table.qza

fi

# Export the relative feature sequence summary table
if [[ ! -d $QIIME_DIR/${ID}_qiime_relative_feature_table ]]; then

qiime tools export \
--input-path $QIIME_DIR/${ID}_relative_feature_table.qza \
--output-path $QIIME_DIR/${ID}_qiime_relative_feature_table

fi

# ------------------------------------------
# ------------------------------------------
# Assign taxonomies to each feature
if [[ ! -f $QIIME_DIR/${ID}_taxonomy.qza ]]; then

qiime feature-classifier classify-consensus-blast \
--i-query $QIIME_DIR/${ID}_feature_sequences.qza \
--i-reference-reads $SILVA_SEQ_DB \
--i-reference-taxonomy $SILVA_TAXONOMY_DB \
--o-classification $QIIME_DIR/${ID}_taxonomy.qza \
--verbose

fi

# Export the taxonomy assignments
if [[ ! -d $QIIME_DIR/${ID}_qiime_taxonomy ]]; then
qiime tools export \
--input-path $QIIME_DIR/${ID}_taxonomy.qza \
--output-path $QIIME_DIR/${ID}_qiime_taxonomy
fi
# ------------------------------------------

# ------------------------------------------
# Export the barplot visualisation artifact
if [[ ! -d $QIIME_DIR/${ID}_qiime_barplot ]]; then

qiime taxa barplot \
--i-table $QIIME_DIR/${ID}_feature_table.qza \
--i-taxonomy $QIIME_DIR/${ID}_taxonomy.qza \
--m-metadata-file $QIIME_DIR/${ID}_sample_metadata.tsv \
--o-visualization $QIIME_DIR/${ID}_taxonomy_barplot.qzv

qiime tools export \
--input-path $QIIME_DIR/${ID}_taxonomy_barplot.qzv \
--output-path  $QIIME_DIR/${ID}_qiime_barplot

fi
# ------------------------------------------

# ------------------------------------------------------------
# Now create the final count and abundance tables

if [[ ! -f $QIIME_DIR/feature_statistics.csv ]]; then
# First transpose the feature tables
if [[ ! -f $QIIME_DIR/${ID}_feature_table-transposed.qza ]]; then
qiime feature-table transpose \
    --i-table $QIIME_DIR/${ID}_feature_table.qza \
    --o-transposed-feature-table $QIIME_DIR/${ID}_feature_table-transposed.qza
fi

# Then combine the transposed feature table with the taxonomy and feature sequences
qiime metadata tabulate  \
--m-input-file $QIIME_DIR/${ID}_taxonomy.qza \
--m-input-file $QIIME_DIR/${ID}_feature_table-transposed.qza \
--m-input-file $QIIME_DIR/${ID}_feature_sequences.qza \
--output-dir $QIIME_DIR/temp

# Export the combined tables
qiime tools export --input-path $QIIME_DIR/temp/visualization.qzv --output-path $QIIME_DIR/temp/out

# Clean up tables
head -n 1 $QIIME_DIR/temp/out/metadata.tsv | sed "s/id\t/ASV\t/g" | sed "s/Consensus/Confidence/g" | sed "s/\t/,/g" > $QIIME_DIR/feature_statistics.csv
tail -n +3 $QIIME_DIR/temp/out/metadata.tsv | sed "s/\t/,/g" >> $QIIME_DIR/feature_statistics.csv

# Remove temporary directories and files
rm -r $QIIME_DIR/temp
rm $QIIME_DIR/${ID}_feature_table-transposed.qza
fi
# ------------------------------------------------------------

# ------------------------------------------------------------
# Now do the same for the relative abundance table. We need to do this a bit differently
# because QIIME can't tranpose a relative frequency table
if [[ ! -f $QIIME_DIR/relative_feature_statistics.csv ]]; then
echo -e "#OTUID\ttaxonomy\tconfidence" > $QIIME_DIR/${ID}_qiime_taxonomy/taxonomy_modified.tsv
tail -n +2 $QIIME_DIR/${ID}_qiime_taxonomy/taxonomy.tsv >> $QIIME_DIR/${ID}_qiime_taxonomy/taxonomy_modified.tsv

biom add-metadata -i $QIIME_DIR/${ID}_qiime_relative_feature_table/feature-table.biom \
-o $QIIME_DIR/${ID}_relative_feature-table_with_taxonomy.biom \
--observation-metadata-fp $QIIME_DIR/${ID}_qiime_taxonomy/taxonomy_modified.tsv \
--sc-separated taxonomy

biom convert -i $QIIME_DIR/${ID}_relative_feature-table_with_taxonomy.biom \
-o $QIIME_DIR/${ID}_relative_feature-table_with_taxonomy.tsv --to-tsv --header-key taxonomy

qiime metadata tabulate \
--m-input-file $QIIME_DIR/${ID}_qiime_taxonomy/taxonomy.tsv \
$QIIME_DIR/${ID}_relative_feature-table_with_taxonomy.tsv \
--m-input-file $QIIME_DIR/${ID}_feature_sequences.qza \
--output-dir $QIIME_DIR/temp2

qiime tools export \
--input-path $QIIME_DIR/temp2/visualization.qzv \
--output-path $QIIME_DIR/temp2/out

head -n 1 $QIIME_DIR/temp2/out/metadata.tsv | sed "s/id\t/ASV\t/g" | sed "s/Consensus/Confidence/g" | sed "s/\t/,/g" > $QIIME_DIR/relative_feature_statistics.csv
tail -n +3 $QIIME_DIR/temp2/out/metadata.tsv | sed "s/\t/,/g" >> $QIIME_DIR/relative_feature_statistics.csv

# Remove duplicate taxonomy column
taxonomy_column=$(head -n 1 $QIIME_DIR/relative_feature_statistics.csv | sed "s/,/\n/g" | grep -n "taxonomy" | cut -d ":" -f 1)
awk -F "," -v tax_col=$taxonomy_column 'BEGIN{OFS=","}{$tax_col=""; gsub(FS "+",FS); print}' $QIIME_DIR/relative_feature_statistics.csv > $QIIME_DIR/temp
mv $QIIME_DIR/temp $QIIME_DIR/relative_feature_statistics.csv

# Remove temporary directories and files
rm -r $QIIME_DIR/temp2
rm $QIIME_DIR/${ID}_relative_feature-table_with_taxonomy.biom
rm $QIIME_DIR/${ID}_relative_feature-table_with_taxonomy.tsv

fi
# ------------------------------------------------------------

exit 1
# HAVE NOT TESTED below

qiime alignment mafft --i-sequences $QIIME_DIR/${ID}_feature_sequences.qza --o-alignment $QIIME_DIR/${ID}_aligned_feature_sequences.qza
qiime alignment mask --i-alignment $QIIME_DIR/${ID}_aligned_feature_sequences.qza --o-masked-alignment $QIIME_DIR/${ID}_masked_aligned_feature_sequences.qza
qiime phylogeny fasttree --i-alignment $QIIME_DIR/${ID}_masked_aligned_feature_sequences.qza --o-tree $QIIME_DIR/${ID}_unrooted_tree.qza
qiime phylogeny midpoint-root --i-tree $QIIME_DIR/${ID}_unrooted_tree.qza --o-rooted-tree $QIIME_DIR/${ID}_rooted_tree.qza

qiime diversity core-metrics --i-phylogeny $QIIME_DIR/${ID}_rooted_tree.qza \
--i-table $QIIME_DIR/${ID}_feature_table.qza \
--p-sampling-depth 5000 \
--output-dir $QIIME_DIR/${ID}_core_metrics

exit 1


# Summarise the features
qiime feature-table summarize \
--i-table $QIIME_DIR/${ID}_feature_table.qza \
--o-visualization $QIIME_DIR/${ID}_feature_table.qzv

# Export the feature table summary
qiime tools export \
--input-path $QIIME_DIR/${ID}_feature_table.qzv \
--output-path $QIIME_DIR/${ID}_qiime_feature_table

# Output the visualisation artifact for the feature sequences
qiime feature-table tabulate-seqs \
--i-data $QIIME_DIR/${ID}_feature_sequences.qza \
--o-visualization $QIIME_DIR/${ID}_feature_seqs.qzv