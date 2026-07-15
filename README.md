# TRITEX long-read assembly pipeline

Table of Contents

* Introduction
  + Installing the required tools
  + Setting up the folder structure
* 1. Obtaining the input contig assembly
* 2. Mapping Hi-C data
* 3. Creating the guide map table
  + 3.1. If you already have a guide map
  + 3.2. If you are using a reference genome to create the guide map
  + 3.3. Mapping the guide map sequences to the genome
* 4. Creating the assembly object
  + 4.1. Breaking chimeric scaffolds
* 5. Manual curation of scaffolds
  + 5.1. Collinearity plots
  + 5.2. Checking Hi-C maps
* 6. Compiling pseudomolecules

This is the documentation of TRITEX, a computational pipeline for chromosome-scale assembly of plant genomes. It was developed in the research group Domestication Genomics at the Leibniz Institute of Plant Genetics and Crop Research (IPK) Gatersleben.

<img width="431" height="135" alt="logo" src="https://github.com/user-attachments/assets/5d95e4bd-5261-4d8b-b4eb-3c8b14bd3401" />

[Legal notice](https://www.ipk-gatersleben.de/en/imprint)

## Introduction

The purpose of the TRITEX long-read assembly pipeline is to build chromosome-scale sequence assemblies for inbred plant genomes. As such, it also applicable in pan-genome projects.

TRITEX uses an input contig assembly (using [HiFi long-reads](https://www.nature.com/articles/s41587-019-0217-9), for example), Hi-C reads and a guide map to assemble the genomes within a short time frame. Manual curation of contig placements is done intuitively with user-editable tables and plots.


The TRITEX workflow can be broadly divided into two stages:

1. *Steps 1-3 of this tutorial.* You will run shell scripts combining standard bioinformatics tools into pipelines for processing Hi-C reads and aligning guide map markers.
2. *Steps 4-6 of this tutorial.* Outputs of phase 1 are loaded into R to interactively create a TRITEX assembly object with tables listing Hi-C links and guide map alignment records.


<img width="2220" height="500" alt="tritex_overview" src="https://github.com/user-attachments/assets/aa5abc66-4efa-4e2e-962e-4869773ed3a0" />
Figure 1. Graphical overview of the TRITEX pipeline.

The core algorithm for Hi-C map construction searches for a minimum spanning tree in the graph induced by Hi-C contact matrix and further refines it to include as many scaffolds as possible and to orient them relative to the chromosomal orientation of the guide map. The algorithm has been described in detail by [Beier *et al.* (2017)](https://doi.org/10.1038/sdata.2017.44:).

In case of questions on the pipeline, contact [Martin Mascher](mailto:mascher@ipk-gatersleben.de).


### Installing the required tools

1. The TRITEX pipeline needs to be run on rather powerful Unix servers. It has been tested on servers running recent CentOS Linux version controlled by the SLURM resource scheduler. Some commands like Hi-C read mapping can be submmitted as batch job. Pseudomolecule construction in R requires interactive sessions, either by direct SSH access or use of interactives SLURM sessions. In any case, the use of a terminal multiplexer such as [tmux](https://github.com/tmux/tmux/wiki) is recommended.
   Some jobs may run for several days. If [Kerberos](https://web.mit.edu/kerberos/) authentification is used, make sure that tickets are renewed.
2. Install the following software. The pipeline has been run successfully with the versions indicated below. Using the most recent versions should be fine as well.

   1. [Z shell](http://zsh.sourceforge.net/), 5.0.2
   2. [GNU Parallel](https://www.gnu.org/software/parallel/), 20150222, [paper](https://www.usenix.org/publications/login/february-2011-volume-36-number-1/gnu-parallel-command-line-power-tool)
   3. [Hifiasm](https://hifiasm.readthedocs.io), 0.15.1-r334, [paper](https://www.nature.com/articles/s41592-020-01056-5)
   4. [gfatools](https://https://github.com/lh3/gfatools), 0.5
   5. [SAMtools](https://github.com/samtools/samtools), 1.9, [paper](https://academic.oup.com/bioinformatics/article/25/16/2078/204688)
   6. bgzip and tabix, part of [BCFtools](https://github.com/samtools/bcftools), 1.9, [paper](https://academic.oup.com/bioinformatics/article/27/5/718/262743)
   7. [Seqtk](https://github.com/lh3/seqtk), 1.0-r76
   8. [bedtools](https://github.com/arq5x/bedtools2), 2.26.0, [paper](https://academic.oup.com/bioinformatics/article/26/6/841/244688)
   9. [R](https://www.r-project.org/), 3.5.1, [paper](https://amstat.tandfonline.com/doi/abs/10.1080/10618600.1996.10474713#.XDDdTy2ZP1I)
   10. [pigz](https://zlib.net/pigz/), 2.3.4, parallel gzip
   11. [minimap2](https://github.com/lh3/minimap2), 2.14, [paper](https://academic.oup.com/bioinformatics/article-abstract/34/18/3094/4994778?redirectedFrom=PDF)
   12. [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html), 1.15, [paper](http://journal.embnet.org/index.php/embnetjournal/article/view/200)
   13. [EMBOSS](http://emboss.sourceforge.net/), 6.6.0, [paper](https://doi.org/10.1016/S0168-9525(00)02024-2)
   14. [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/), 37.93
   15. faToTwoBit and twoBitInfo, from the [UCSC Genome Browser](http://hgdownload.cse.ucsc.edu/admin/exe/), v.3.69
   16. [Novosort](http://www.novocraft.com/products/novosort/), 3.06.05. This is commercial software. A single-threaded version is available for free. Novosort is used only in Hi-C read mapping scripts. An alternative using SAMtools instead is available [here](https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/shell/).

   The pipeline assumes that the GNU (as opposed to BSD) versions of common UNIX tools (`coreutils`, `sed`, `awk`) are used.
3. The following R packages should be installed and accessible in `.libPaths()`. The versions specified are known to work with TRITEX. Newer versions will likely work as well.

   1. [data.table](https://cran.r-project.org/web/packages/data.table/index.html), 1.11.8
   2. [stringi](https://cran.r-project.org/web/packages/stringi/index.html), 1.2.4
   3. [igraph](https://cran.r-project.org/web/packages/igraph/index.html), 1.2.2, [paper](https://static1.squarespace.com/static/5b68a4e4a2772c2a206180a1/t/5cd1e3cbb208fc26c99de080/1557259212150/c1602a3c126ba822d0bc4293371c.pdf)
   4. [zoo](https://cran.r-project.org/web/packages/zoo/index.html), 1.8-3, [paper](https://arxiv.org/abs/math/0505527)
   5. [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf), included in R-core in R versions > 2.14
4. If you want to use the [R Shiny](https://shiny.rstudio.com/) app Hi-C map inspector, you need to set up a web server running R Shiny. Refer to the R Shiny [manual](https://www.linode.com/docs/guides/how-to-deploy-rshiny-server-on-ubuntu-and-debian/). Note that setting up a Shiny server is a non-trivial task of Unix systems administration and requires super-user rights.
5. Here are some thoughts on hardware requirements and how to allocate resource on institutional resource scheduling systems:

   1. For the hardware requirements of primary contig assembly, we refer to the documentation of hifiasm[<https://github.com/chhylp123/hifiasm#results>], a tool commonly used for that purpose. Thirty 30 and 1 TB RAM will suffice to assemble genomes sequences of 15 GB size within a few days.
   2. The Hi-C read mapping script uses `minimap2`, which is multi-threaded.
      More cores translate into shorter run times approximately linearly.
      Sorting of alignment records `Novosort` (and `samtools`) can be speeded up by allocating more RAM for in-memory sorting. Otherwise, temporary disk space is used (which incurs longer run-time because of disk I/O).
   3. The interactive R sessions for pseudomolecule construction hold large tables of Hi-C links in main memory. For a genome in size range of 5-15 Gb, 500 GB - 1 TB of main memory may be required.
      Some parts of the pseudomolecule construction are multithreaded (e.g. the `hic_map()` function),
      but savings in run time are negligeable, so running R on an allocation of 1 core will be just fine.
6. Clone the Bitbucket [repository](https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/) of the TRITEX assembly pipeline. The path to the local copy of the repository will be referenced by `$bitbucket` in the code listings.
7. The scripts used below have parameters such as `--bedtools` or `--samtools` (indicated in the code listings) to specify the absolute paths to the required executables. For better readability, they are omitted from the commands.


| Instead of specifying the paths each time you call the scripts, you may also edit the scripts and modify the default paths. We chose not to use executables in `$PATH` to enforce proper documentation of software versions. |

### Setting up the folder structure

A folder for each assembly project should be created. We refer to this directory as `$projectdir` in the code listings below.

The project folder should have the following subdirectory structure. At the start of the pipeline, all folders except "raw\_reads" should be empty.

* **raw\_reads**

  + **hifi\_reads**
  + **hic\_reads**

    - These folders contain the raw reads (or symbolic links to them).
* **assembly**

  + This folder will contain the CCS assembly.
* **mapping**

  + This folder will contain the results of mapping of Hi-C reads to the HiFi assembly.

    - **hic**

      * In this folder you should create symbolic links to the Hi-C reads.
* **pseudomolecules**

  + This folder will be used to run the Pseudomolecule construction section of this tutorial. It will contain the R objects of the assemblies, plots for quality assessment and final pseudomolecule FASTA files.

## 1. Obtaining the input contig assembly

1. Any long-read input assembly will work, but in this tutorial we are going to use HiFi reads with hifiasm.
2. Let’s assemble the HiFi reads with hifiasm. Set up a list of all FASTQ files and run hifiasm. Results will be inside the `$projectdir/assembly/project_name/` directory.

   ```
   cd $projectdir/assembly
   d='$projectdir/raw_reads/hifi_reads'
   find $d | grep 'fastq.gz$' > hifi_reads.txt

   reads='hifi_reads.txt'
   prefix='project_name'
   threads=30 (1)
   
   mkdir $prefix
   out=$prefix/$prefix

   xargs -a $reads hifiasm -t $threads -o $out > $out.out 2> $out.err
   ```

   **1** | Set the number of threads you would like to use. |
3. Check the stats of contigs and unitigs.

   ```
   cd $projectdir/assembly/project_name/

   grep '^S' project_name.p_ctg.noseq.gfa | cut -f 4 | cut -d : -f 1,3 \
    | tr : '\t' | $bitbucket/shell/n50 > project_name.p_ctg.noseq.gfa.n50 (1)

   grep '^S' project_name.p_utg.noseq.gfa | cut -f 4 | cut -d : -f 1,3 \
   | tr : '\t' | $bitbucket/shell/n50 > project_name.p_utg.noseq.n50
   ```

   **1** | The backslash \ will be used throughout the tutorial to break long shell commands into several lines. When pasting the commands or editing them, make sure that no white space follows the backslash. Otherwise, the shell will interpret the lines as belonging to different commands. Also multi-line commands do not tolerate intervening command line (starting with the hash sign #). |
4. Convert contig GFA to FASTA and change contig names.

   ```
   gfa='project_name.p_ctg.gfa'
   gfatools gfa2fa $gfa | seqtk rename - contig_ > ${gfa:r}.fa (1)
   ```

   **1** | The syntax `${gfa:r}` is specific to the Z shell, see the section of the shell manual on [expansion and substitution](https://linux.die.net/man/1/zshexpn). Very usefyl are `:r` and `:t` modifiers to remove filename extensions and leading pathname components, respectively. The latter works like the UNIX command `basename`. |

## 2. Mapping Hi-C data

1. THe first step is to digest *in silico* the hifiasm assembly with the restriction enzyme used for the preparation of the Hi-C libraries.

   ```
   cd $projectdir/assembly/project_name/
   ref='$projectdir/assembly/project_name/project_name.p_ctg.fa' (1)

   $bitbucket/shell/digest_emboss.zsh --ref $ref --enzyme 'MboI' --sitelen 4 --minlen 30  (2) (3)
   ```

   **1** | `$ref` is the renamed FASTA created in the previous step. |
   **2** | In the script "digest\_emboss.zsh", change the paths to the executables of `restrict`, `bedtools` and `rebase`. |
   **3** | Don’t forget to make sure that the enzyme matches the one used for making your Hi-C libraries. Another common choice is `--enzyme 'DpnII' --sitelen 4`. |
2. Put the gzipped FASTQ files of the Hi-C reads (or symbolic links to them) of the Hi-C reads in the folder "hic".
3. Now start the Hi-C mapping.

```
cd $projectdir/mapping/
ref='$projectdir/assembly/project_name/project_name.p_ctg.fa'
map='$bitbucket/shell/run_hic_mapping.zsh' (1)
bed=${ref:r}_MboI_fragments_30bp.bed (2)
cd $projectdir

$map --threads 64 --mem '200G' --linker "GATCGATC" --ref $ref --bed $bed --tmp $TMPDIR hic (3)(4)
```

**1** | You need to specify in the script the paths to the following executables: Cutadapt (`--cutadapt`), bgzip (`--bgzip`), Minimap2 (`--minimap`), Novosort (`--novosort`), SAMtools (`--samtools`), and BEDTools (`--bedtools`). |
**2** | Use the BED file with the in silico digest for the correct restriction enzyme. Using the wrong enzyme will give wrong results. |
**3** | Specify correct linker sequences that is created by ligation of fill-in overhangs during the Hi-C protocol. If HindIII was used, the correct linker is AAGCTAGCTT; if DpnII was used, the linker is GATCGATC. |
**4** | The temporary directory needs to be large enough to hold huge BAM files (dozens of GB). |

## 3. Creating the guide map table

### 3.1. If you already have a guide map

1. Follow these steps if you already have a guide map. We are going to adjust it to the right format to be used in TRITEX.
2. If you plan on using a reference genome to build a guide map, go to the next section.
3. An example file is available [here](https://doi.org/10.5447/ipk/2022/20) (MaizeB73\_guidemap\_marker\_IBM.Rds).
4. You will need a BED file with the genomic positions of your markers. The fourth column should contain IDs like "seq\_N", where N goes from 1 to the number of markers.
5. Most steps will be done in R from now on. If you are unsure whether a command should be pasted to the R or the shell prompt, hover with the mouse over the code listings. SH or R should appear in the top-right corner.

   ```
   # read the BED file containing the marker regions
   fread('pseudomolecules_singlecopy_100bp.bed') -> p

   setnames(p, c("agp_chr", "agp_start", "agp_end", "seq"))
   p[, .(css_contig=seq, popseq_cM=(agp_start+agp_end)/2/1e6, #set pseudo-cM as Mb position
       sorted_genome = as.character(NA),
       css_contig_length = agp_end - agp_start,
       popseq_alphachr = as.character(agp_chr), # it needs to be character column
       sorted_arm = NA) ] -> pp

   pp[,popseq_chr:=as.integer(popseq_alphachr)]
   pp[, sorted_alphachr := popseq_alphachr]
   pp[, sorted_chr := popseq_chr]
   pp[, sorted_lib := popseq_alphachr]

   # Save as RDS file for use in TRITEX pipeline
   saveRDS(pp, "project_name_pseudopopseq.Rds")
   ```

### 3.2. If you are using a reference genome to create the guide map

1. If you have a guide map available, you can skip this step and go to the next section.
2. In case you want to create a guide map from an available reference genome, first you need to obtain single-copy 100bp regions.

```
mask='$bitbucket/miscellaneous/mask_assembly.zsh' (1)
fa='$projectdir/reference_assembly.fasta' (2)
bed="${fa:t:r}_masked_noGaps.bed"
out="${fa:t:r}_singlecopy_100bp"

$mask --mem 300G --fasta $fa --mincount 2 --out "."

awk '$3 - $2 >= 100' $bed | grep -v chrUn | awk '{print $0"\tseq_"NR}' \
 | tee $out.bed | bedtools getfasta -fi $fa -bed /dev/stdin -name -fo $out.fasta
```

**1** | You need to specify the paths to the following executables: BBDuk (`bbduk`), fatotwobit (`bit`), twobitinfo (`tbinfo`), kmercountexact.sh (from BBMap) (`kmer`), SAMtools (`samtools`), and BEDTools (`bedtools`). |
**2** | Use a reference assembly for the studied species. |

1. Most steps will be done in R from now on. If you are unsure whether a command should be pasted to the R or the shell prompt, hover with the mouse over the code listings. SH or R should appear in the top-right corner.
2. Now it is time to build the pseudo POPSEQ file (the guide map), which is used as input for the TRITEX pipeline.
3. [Here](https://doi.org/10.5447/ipk/2022/20) is an example file of how your guide map should look like.

```
# Read single-copy regions generated in last step
f<-'reference_assembly_singlecopy_100bp.bed'
fread(f) -> p
setnames(p, c("agp_chr", "agp_start", "agp_end", "seq"))

# Change chromosome names to their respective chr numbers
# These are example names from maize

ids <- c("NC_050096.1",
         "NC_050097.1",
         "NC_050098.1",
         "NC_050099.1",
         "NC_050100.1",
         "NC_050101.1",
         "NC_050102.1",
         "NC_050103.1",
         "NC_050104.1",
         "NC_050105.1")

for (i in 1:10){
  p[grepl(ids[i], agp_chr), agp_chr:=as.integer(i)]
}

# This is for unplaced scaffolds, assigning them as "NA"
p[grepl("NW", agp_chr), agp_chr:=as.integer(NA)]

p[, .(css_contig=seq, popseq_cM=(agp_start+agp_end)/2/1e6, #set pseudo-cM as Mb position
    sorted_genome = as.character(NA),
    css_contig_length = agp_end - agp_start,
    popseq_alphachr = as.character(agp_chr), # it needs to be character column
    sorted_arm = NA) ] -> pp

pp[,popseq_chr:=as.integer(popseq_alphachr)]
pp[, sorted_alphachr := popseq_alphachr]
pp[, sorted_chr := popseq_chr]
pp[, sorted_lib := popseq_alphachr]

# Save as RDS file for use in TRITEX pipeline
saveRDS(pp, "project_name_pseudopopseq.Rds")
```

### 3.3. Mapping the guide map sequences to the genome

1. Independently of how you obtained your guide map, you should map the marker sequences to the contig assembly.
2. The output PAF file will be input into R into the next step.

```
minimap='/opt/Bio/minimap2/2.17/bin/minimap2'
ref='$projectdir/assembly/project_name/project_name.p_ctg.fa'
qry='reference_assembly_singlecopy_100bp.fasta' # this is the guide map sequence file
size='20G'
mem='5G'
threads=30
prefix=${ref:t:r}_project_name

{
$minimap -t $threads -2 -I $size -K$mem -x asm5 $ref $qry | bgzip > $prefix.paf.gz
} 2> $prefix.err &
```

## 4. Creating the assembly object

1. Now we are going to load the contig assembly, guide map and the Hi-C alignment records into R to create the assembly object.
2. All subsequent steps should be carried out in the working directory `$projectdir/pseudomolecules/`.
3. First load all the necessary datasets into R.

   ```
   # Load the R functions for pseudomolecule construction.
   source("$bitbucket/R/pseudomolecule_construction.R") (1)

   # Import the guide map.
   readRDS('$projectdir/project_name_pseudopopseq.Rds') -> popseq

   # Read contig lengths
   f <- '$projectdir/assembly/project_name/project_name.p_ctg.fa.fai'
   fread(f, head=F, select=1:2, col.names=c("scaffold", "length")) -> fai

   # Read guide map alignment.
   f <- '$projectdir/assembly/project_name/project_name_ctg.paf.gz' (2)
   read_morexaln_minimap(paf=f, popseq=popseq, minqual=30, minlen=500, prefix=F) -> morexaln (3)

   # Read Hi-C links
   dir <- '$projectdir/hic'
   fread(paste('find', dir, '| grep "fragment_pairs.tsv.gz$" | xargs zcat'),
          header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2")) -> fpairs
   ```

   **1** | Note that `$bitbucket`, `$projectdir` and `$assembly` will not be evaluated as variables by R here and in all other instances below. Replace them always with a correct path by modifying the command with a text editor before running it. |
   **2** | Note that IDs in this file should be the same as in the guide map. Modify them if needed. |
   **3** | Change the `minlen` parameter depending on the length of the guide map marker size. |
4. Now initialize and save the assembly.

   ```
   # Init and save assembly, replace species with the right one. Here we are going to use "maize".
   init_assembly(fai=fai, cssaln=morexaln, fpairs=fpairs) -> assembly
   anchor_scaffolds(assembly = assembly, popseq=popseq, species="maize") -> assembly (1)
   add_hic_cov(assembly, binsize=1e4, binsize2=1e6,
   	minNbin=50, innerDist=3e5, cores=40) -> assembly  (2) (3)
   saveRDS(assembly, file="assembly.Rds") (4)
   ```

   **1** | Possible values for the species parameter are: 'barley', 'wheat' (works for A genome species, *T. turgidum* and *T. aestivum*), 'maize', 'rye', 'oats\_new' (chromosomes 1 to 7, subgenomes A, B and D), 'avena\_barbata', 'hordeum\_bulbosum', 'faba\_bean', and 'lolium'. More species can be accomodated by re-defining the function `chrNames()`. Type chrNames to see its source code. The only use of species is to map numeric chromosome IDs to proper names, e.g. 1 to 1H in barley (see the note by [Linde-Laursen](https://wheat.pw.usda.gov/ggpages/bgn/26/text261a.html) to see why this makes a difference) or 21 to 7D. |
   **2** | Change number of cores according to match your allocated resources. Don’t forget to change it in the next steps also. |
   **3** | `binsize` is the resolution at which physical coverage is calculated. `binsize2` is the size of batches for parallelization; smaller values decrease memory consumption. |
   **4** | It is recommended to prefix the names of the RDS files created in this and later steps with the current date (e.g. "220801\_assembly.Rds") to distinguish different versions. |

### 4.1. Breaking chimeric scaffolds

1. Generate diagnostic plots to check for chimeras in the input assembly. If sequences from different chromosomes are joined, the chimeric scaffolds will have guide map markers from more than one chromosome.

   ```
   # Create diagnostic plots for contigs longer than 1 Mb.
   assembly$info[length >= 1e6, .(scaffold, length)][order(-length)] -> s
   plot_chimeras(assembly=assembly, scaffolds=s, species="maize",
   	refname="B73", autobreaks=F, mbscale=1, file="assembly_1Mb.pdf", cores=40) (1)

   # Diagnostic plots for putative chimeras with strong drops in Hi-C coverage.
   assembly$info[length >= 1e6 & mri <= -3, .(scaffold, length)][order(-length)] -> s
   plot_chimeras(assembly=assembly, scaffolds=s, species="maize",
   	refname="B73", autobreaks=F, mbscale=1, file="assembly_chimeras.pdf", cores=40)
   ```

   **1** | Don’t forget to change species and reference name. |
2. Check the plots and search for drops in Hi-C coverage, misalignments and chimeras. For each case, assign the PDF page and scaffold position which needs breaking.
3. For instance, in [this example](https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/raw/9375957ff5f1763b1ce11d090919a76de9d7bf7a/example_assembly_1Mb.pdf), the scaffolds that need to be broken are on pages 3, 6, and 9. We need to change the bin region on the code for each scaffold, like this:

   ```
   assembly$info[length >= 1e6, .(scaffold, length)][order(-length)] -> ss

    i=3 # page number in PDF
    ss[i]$scaffold -> s
    # change the bin position observed in the plots
    assembly$cov[s, on='scaffold'][bin >= 10e6 & bin <= 22e6][order(r)][1, .(scaffold, bin)] -> b

    i=6 # page number in PDF
    ss[i]$scaffold -> s
    rbind(b, assembly$cov[s, on='scaffold'][bin >= 40e6 & bin <= 45e6][order(r)][1, .(scaffold, bin)])->b

    i=9 # page number in PDF
    ss[i]$scaffold -> s
    rbind(b, assembly$cov[s, on='scaffold'][bin >= 1e6 & bin <= 9e6][order(r)][1, .(scaffold, bin)])->b

   # And so on.

   # Rename column and plot again to double-check the breakpoints are correct
   setnames(b, "bin", "br")
   plot_chimeras(assembly=assembly, scaffolds=b, br=b, species="maize",
   	refname="B73",  mbscale=1, file="assembly_chimeras_final.pdf", cores=40)
   ```
4. After checking if the break points are in the right places, proceed with breaking the scaffolds.

   ```
   break_scaffolds(b, assembly, prefix="contig_corrected_v1_", slop=1e4, cores=40, species="maize") -> assembly_v2
   saveRDS(assembly_v2, file="assembly_v2.Rds") (1)
   ```

   **1** | This will create a new file for the assembly. |
5. Now we make the Hi-C map. Three files will be generated: a PDF with the intrachromosomal plots, a PDF with the interchromosomal plots, and a file to be loaded into R Shiny Map inspector (`*_export.Rds`).

   ```
   fbed <- '$projectdir/assembly/project_name/project_name.p_ctg_MboI_fragments_30bp.bed'
   read_fragdata(info=assembly_v2$info, file=fbed) -> frag_data

   # Consider only contigs >= 1 Mb first
   frag_data$info[!is.na(hic_chr) & length >= 1e6, .(scaffold, nfrag, length, chr=hic_chr, cM=popseq_cM)] -> hic_info

   hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed,
   	species="maize", ncores=40, min_nfrag_scaffold=30,
   	max_cM_dist = 1000, binsize=2e5, min_nfrag_bin=10, gap_size=100) -> hic_map_v1

   saveRDS(hic_map_v1, file="hic_map_v1.Rds")

   # Make the Hi-C plots, files will be placed in your working directory.
   snuc <- '$projectdir/assembly/project_name/project_name.p_ctg_MboI_fragments_30bp_split.nuc.txt'

   hic_plots(rds="hic_map_v1.Rds", assembly=assembly_v2,
   	cores=40, species="maize", nuc=snuc) -> hic_map_v2
   ```
6. Let’s do a second round of checking for chimeras. Broken scaffolds in the last step are renamed in the plots. It might not be necessary to break any more scaffolds; in this case, just check the plots and if everything is alright, proceed to step 8.

   ```
   # If you'd like to  plot specific contigs:
   data.table(scaffold=c('contig_corrected_v1_53', 'contig_corrected_v1_24',
   	'contig_corrected_v1_5', 'contig_corrected_v1_76') -> ss

   # If you'd like to check them all again:
   assembly_v2$info[length >= 1e6, .(scaffold, length)][order(-length)] -> ss

   plot_chimeras(assembly=assembly_v2, scaffolds=ss,
           species="maize", refname="B73", autobreaks=F,
           mbscale=1, file="assembly_v2_chimeras.pdf", cores=40)

   # In case there are more chimeras, proceed like the previous step. Assign `i` with the PDF page and don't forget to change the bin.

   i=1 # PDF page
   ss[i]$scaffold -> s
   assembly_v2$cov[s, on='scaffold'][bin >= 6e7 & bin <= 8e7][order(r)][1, .(scaffold, bin)] -> b

   i=2 # PDF page
   ss[i]$scaffold -> s
   rbind(b, assembly_v2$cov[s, on='scaffold'][bin >= 8e7 & bin <= 10e7][order(r)][1, .(scaffold, bin)]) -> b

   setnames(b, "bin", "br")
   plot_chimeras(assembly=assembly_v2, scaffolds=b,
   	br=b, species="maize", refname="B73",  mbscale=1,
   	file="assembly_v2_chimeras_final.pdf", cores=40)

   break_scaffolds(b, assembly_v2, prefix="contig_corrected_v2_", slop=1e4,
   	cores=40, species="maize") -> assembly_v3 (1)
   saveRDS(assembly_v3, file="assembly_v3.Rds") (1)

   fbed <- '$projectdir/assembly/project_name/project_name.p_ctg_MboI_fragments_30bp.bed'
   read_fragdata(info=assembly_v3$info, file=fbed) -> frag_data

   frag_data$info[!is.na(hic_chr) & length >= 1e6, .(scaffold, nfrag, length, chr=hic_chr, cM=popseq_cM)] -> hic_info

   hic_map(info=hic_info, assembly=assembly_v3, frags=frag_data$bed,
   	species="maize", ncores=40, min_nfrag_scaffold=30, max_cM_dist = 1000,
   	binsize=2e5, min_nfrag_bin=10, gap_size=100) -> hic_map_v2 (1)

   saveRDS(hic_map_v2, file="hic_map_v2.Rds") (1)

   snuc <- '$projectdir/assembly/project_name.p_ctg_MboI_fragments_30bp_split.nuc.txt'
   hic_plots(rds="hic_map_v2.Rds", assembly=assembly_v3,
   	cores=40, species="maize", nuc=snuc) -> hic_map_v2 (1)
   ```

   **1** | Don’t forget to update assembly and Hi-C map versions. |
7. If another round of breaking scaffolds is needed, you can just run the previous block again. Don’t forget to change the output files and object names (assembly\_vX, hic\_map\_vX).
8. If no more chimeras are found, proceed to decreasing the size to 500 kb.

   ```
   frag_data$info[!is.na(hic_chr) & length >= 5e5, .(scaffold, nfrag, length, chr=hic_chr, cM=popseq_cM)] -> hic_info

   hic_map(info=hic_info, assembly=assembly_v3, frags=frag_data$bed, species="maize",
   	ncores=40, min_nfrag_scaffold=30, max_cM_dist = 1000, binsize=2e5,
   	min_nfrag_bin=10, gap_size=100) -> hic_map_v3 (1)

   saveRDS(hic_map_v3, file="hic_map_v3.Rds") (1)

   snuc <- '$projectdir/assembly/project_name.p_ctg_MboI_fragments_30bp_split.nuc.txt'
   hic_plots("hic_map_v2.Rds", assembly=assembly_v3, cores=40,
   	species="maize", nuc=snuc) -> hic_map_v3 (1)


   # Exporting to Excel the table that will be used for manual curation.
   write_hic_map(rds="hic_map_v3.Rds", file="hic_map_v3.xlsx", species="maize")
   ```

   **1** | Don’t forget to put the correct object name here. |

## 5. Manual curation of scaffolds

Now that we have the scaffolded assembly, we are going to check if there are any inversions, extra contigs, or misplaced contigs. To do so, we can check the Hi-C plots generated in the last step coupled with a collinearity plot of scaffolds (the AGP, `hic_map` object) to the guide map. It is also possible to use the R Shiny Map Inspector app to zoom in the Hi-C intrachromosomal plots to check the contigs that need to be changed.

### 5.1. Collinearity plots

1. This step is optional, yet recommended. It will help spot inverted or misplaced contigs. The collinearity plots will show marker positions and against the assembly. It is possible to check the alignments and spot inverted/misplaced contigs. You can do this anytime you want to check the collinearity.

   ```
   assembly_v3$cssaln[!is.na(mb), .(scaffold, pos, guide_chr=paste0("chr", "", popseq_alphachr), guide_pos=mb*1e6)] -> x (1)
   hic_map_v3$agp[, .(scaffold, agp_start, agp_end, orientation, agp_chr)][x, on="scaffold"] -> x (1)

   # Lift alignment positions to pseudomolecule positions
   x[orientation == 1, agp_pos :=  agp_start + pos]
   x[orientation == -1, agp_pos :=  agp_end - pos]

   # Chromosome-wise plot
   x[agp_chr == guide_chr] -> x
   pdf("hic_map_v3_vs_assembly_v3.pdf", height=3*3, width=3*3)
   par(mar=c(5,5,3,3))
   lapply(chrNames(species="maize", agp=T)$agp_chr[1:21], function(i){
     print(i)
     x[i, plot(pch=".", main=i, agp_pos/1e6, guide_pos/1e6, xlab="Hi-C AGP (Mb)", ylab="guide (Mb)", las=1, bty='l'), on="agp_chr"]

     # Draw contig boundaries, use correct hic_map version here
     hic_map_v3$agp[agp_chr == i & gap == T][, abline(v=agp_start/1e6, col="gray")] (1)
   })
   dev.off()
   ```

   **1** | Don’t forget to change objects' names. |

### 5.2. Checking Hi-C maps

1. We can use Map Inspector Shiny app or the intrachromosomal plots for manual curation. Genomic regions containing putatively chimeric scaffolds can be clicked on to get their names and pinpoint the sites of culpable misjoins.
2. The source code required for deploying Hi-C Map Inspector on an R Shiny server is contained in `$bitbucket/R/map_inspector.Rmd`. You need to modify the paths to the data directory and the default map object.
3. You’ll want to check for inversions near scaffold borders (grey lines).
4. When modifications are spotted, change the Excel file generated in the last steps. If an inverted contig is spotted, the column "new\_orientation" on the table must be changed (Fig. 2A). On the other hand, if there is a misplaced/extra contig, the row containing it should be move (don’t forget to change the bin number order) (Fig. 2B).

<img width="4267" height="2917" alt="manual_curation" src="https://github.com/user-attachments/assets/0e97ce16-241f-491e-a8f7-adb5c5af0e7e" />
Figure 2. Manual curation in the TRITEX’s correct-map-inspect cycle. (A) The Hi-C contacts show a pattern indicative of an inversion in the terminal contig. The orientation is swapped in the Excel table and a new Hi-C matrix is computed with the updated configuration. The revised Hi-C matrix has fewer off-diagonal signals. (B) Hi-C contacts show a pattern indicative of a misplaced contig. The order of the final two rows is reversed in the Excel table and the Hi-C matrix is computed with the new configuration. The revised Hi-C matrix has fewer off-diagonal signals.

1. After making the changes, import the Excel file again and proceed with mapping.

   ```
   read_hic_map(rds="hic_map_v3.Rds", file="hic_map_v3_edit.xlsx") -> nmap
   diff_hic_map(rds="hic_map_v3.Rds", nmap, species="maize")
   hic_map(species="maize", agp_only=T, map=nmap) -> hic_map_v4

   saveRDS(hic_map_v4, file="hic_map_v4.Rds")

   # Final Hi-C plots. You can check if the inverted scaffolds are correct.
   snuc <- '$projectdir/assembly/project_name.p_ctg_MboI_fragments_30bp_split.nuc.txt'
   hic_plots(rds="hic_map_v4.Rds", assembly=assembly_v3,
   	cores=30, species="maize", nuc=snuc) -> hic_map_v4
   ```
2. If there are still wrong contigs, you can make new modifications to the Excel file and repeat this block until there are no more changes needed.

## 6. Compiling pseudomolecules

1. First let’s get some stats of the pseudomolecules.

   ```
   hic_map_v4$agp[agp_chr != "chrUn" & gap == F][, .("Ncontig" = .N,
   	"N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
   	"max_length (Mb)"=round(max(scaffold_length/1e6),1),
   	"min_length (kb)"=round(min(scaffold_length/1e3),1)), key=agp_chr] -> res
   hic_map_v4$chrlen[, .(agp_chr, "length (Mb)"=round(length/1e6, 1))][res, on="agp_chr"] -> res
   setnames(res, "agp_chr", "chr")

   hic_map_v4$agp[gap == F & agp_chr != "chrUn"][, .("Ncontig" = .N,
   	"N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
   	"max_length (Mb)"=round(max(scaffold_length/1e6),1),
   	"min_length (kb)"=round(min(scaffold_length/1e3),1))][, agp_chr := "1-10"] -> res2
   hic_map_v4$chrlen[, .(agp_chr="1-10", "length (Mb)"=round(sum(length)/1e6, 1))][res2, on="agp_chr"] -> res2
   setnames(res2, "agp_chr", "chr")

   hic_map_v4$agp[gap == F & agp_chr == "chrUn"][, .(chr="un", "Ncontig" = .N,
   	"N50 (Mb)"=round(n50(scaffold_length/1e6), 1),
   	"max_length (Mb)"=round(max(scaffold_length/1e6),1),
   	"min_length (kb)"=round(min(scaffold_length/1e3),1),
   	"length (Mb)"=sum(scaffold_length/1e6))] -> res3

   rbind(res, res2, res3) -> res

   write.xlsx(res, file="hic_map_v4_pseudomolecule_stats.xlsx")
   ```
2. And finally compile the pseudomolecules. The FASTA file will be placed inside the directory assigned in "output".

```
fasta <- '$projectdir/assembly/project_name/project_name.p_ctg.fa'
sink("pseudomolecules_v1.log") (1)
compile_psmol(fasta=fasta, output="pseudomolecules_v1",
	hic_map=hic_map_v4, assembly=assembly_v3, cores=30)
sink() (1)
```

**1** | That opens and closes a log file. |

Last updated 2022-12-15 15:16:29 +0100
