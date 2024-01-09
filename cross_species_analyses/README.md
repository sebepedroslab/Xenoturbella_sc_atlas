# SAMap

Analysis of cross-species transcriptome similarity using SAMap, comparing _Xenoturbella bockii_ (`Xboc`) to four species: _Isodiametra pulchra_ (`Isopu`), _Strongylocentrotus purpuratus_ (`Spur`), _Hofstenia miamia_ (`Hmia`) and _Nematostella vectensis_ (`Nvec`).

Folder structure:

* `data/expression` contains single-cell expression data (per-cell cell type classifications, `mtx` UMI matrices, and gene names).
* `data/reference` contains reference sequence data for each species in the analysis (fasta).
* `data/blast` blast database of pairwise species analyses, formatted for SAMap (gzipped `txt` files).
* `results_alignment_samap`: final results (intermediate SAMap files not included).

Steps:

1. Create UMI matrices (ungzip `mtx` files in `data/expression/` folder as needed):

```R
# Spur, Isopu, Hmia: as is
# Nvec:
library("metacell")
library("Matrix")
metacell::scdb_init("data/expression/", force_reinit = TRUE)
mat = metacell::scdb_mat("nvec")
mc = metacell::scdb_mc("bct_nvec")
sc = data.frame(cell = names(mc@mc), celltype = mc@mc)
sc = sc [ !duplicated(sc$cell), ]
gg = rownames(mat@mat)
gg = gg [ !grepl("orphan_peak", gg) ]
write.table(sc, "data/expression/Nvec_cellID_cell_type", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(gg, "data/expression/Nvec_features.tsv", row.names = FALSE, col.names = FALSE, quote = FALSE)
Matrix::writeMM(mat@mat[gg,sc$cell], "data/expression/Nvec_UMI_table.mtx")
```

2. Create blast databases for SAMap:

```bash
conda activate samap

# blast dbs
makeblastdb -dbtype prot -parse_seqids -in data/reference/Xboc.fasta
makeblastdb -dbtype nucl -parse_seqids -in data/reference/Spur.fasta
makeblastdb -dbtype nucl -parse_seqids -in data/reference/Isopu.fasta
makeblastdb -dbtype nucl -parse_seqids -in data/reference/Hmia.fasta
makeblastdb -dbtype prot -parse_seqids -in data/reference/Nvec.fasta

# blast alignments
# Xboc as reference
bash qsub_blast.sh data/reference/Xboc.fasta data/reference/Spur.fasta  data/samap.Spur-Xboc.blast.csv  blastx
bash qsub_blast.sh data/reference/Xboc.fasta data/reference/Hmia.fasta  data/samap.Hmia-Xboc.blast.csv  blastx
bash qsub_blast.sh data/reference/Xboc.fasta data/reference/Isopu.fasta data/samap.Isopu-Xboc.blast.csv blastx
bash qsub_blast.sh data/reference/Xboc.fasta data/reference/Nvec.fasta  data/samap.Nvec-Xboc.blast.csv  blastp
# Xboc as query
# Xboc as reference
bash qsub_blast.sh data/reference/Spur.fasta  data/reference/Xboc.fasta data/samap.Xboc-Spur.blast.csv  tblsatn
bash qsub_blast.sh data/reference/Hmia.fasta  data/reference/Xboc.fasta data/samap.Xboc-Hmia.blast.csv  tblastn
bash qsub_blast.sh data/reference/Isopu.fasta data/reference/Xboc.fasta data/samap.Xboc-Isopu.blast.csv tblastn
bash qsub_blast.sh data/reference/Nvec.fasta  data/reference/Xboc.fasta data/samap.Xboc-Nvec.blast.csv  blastp

# format blast alignment files for SAMap and wrangle a bit with gene names to fit UMI matrices
# awk 'BEGIN{OFS="\t"}{$1="Spur_"$1; $2="Xboc_"$2 ; print $0}'
mkdir data/blast/XbocSpur
cp data/samap.Spur-Xboc.blast.csv  data/blast/XbocSpur/Spur_to_Xboc.txt    && sed -i "s/\\.i[0-9]*\t/\t/" data/blast/XbocSpur/Spur_to_Xboc.txt
cp data/samap.Xboc-Spur.blast.csv  data/blast/XbocSpur/Xboc_to_Spur.txt    && sed -i "s/\\.i[0-9]*\t/\t/" data/blast/XbocSpur/Xboc_to_Spur.txt
mkdir data/blast/XbocHmia
cp data/samap.Hmia-Xboc.blast.csv  data/blast/XbocHmia/Hmia_to_Xboc.txt
cp data/samap.Xboc-Hmia.blast.csv  data/blast/XbocHmia/Xboc_to_Hmia.txt
mkdir data/blast/XbocIsopu
cp data/samap.Isopu-Xboc.blast.csv data/blast/XbocIsopu/Isopu_to_Xboc.txt  && sed -i "s/_i[0-9]*\t/\t/" data/blast/XbocIsopu/Isopu_to_Xboc.txt
cp data/samap.Xboc-Isopu.blast.csv data/blast/XbocIsopu/Xboc_to_Isopu.txt  && sed -i "s/_i[0-9]*\t/\t/" data/blast/XbocIsopu/Xboc_to_Isopu.txt
mkdir data/blast/XbocNvec
cp data/samap.Nvec-Xboc.blast.csv  data/blast/XbocNvec/Nvec_to_Xboc.txt
cp data/samap.Xboc-Nvec.blast.csv  data/blast/XbocNvec/Xboc_to_Nvec.txt

```
3. Run SAMap for each species pair with `Xboc`:

```bash
python  s31_samap_prepare_SAM_2023-03-20.py
python  s32_samap_objects_2023-03-20.py
python  s33_samap_postalignment_2023-03-20.py
Rscript s34_samap_plots_2023-03-21.R
```






