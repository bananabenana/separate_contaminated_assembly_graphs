# separate_contaminated_assembly_graphs
Script which split assembly graphs into membership groups by taxon. This allows you to decontaminate multi-species assembly graphs into species-specific sub-graphs and assemblies. Will give you seperate `.gfa` and `.fasta` files for each taxa within your contaminated `.gfa` file. You decide the taxonomic level (`"kingdom", "phylum", "class", "order", "family", "genus", "species"`) for which to separate your assemblies.


## Motivation
Assemblies are often chimeric, in that they can be contaminated by multiple species. This is hard to de-convolute using fasta files. Hence, this script works on the connected memberships within the assembly graphs (`.gfa` files). It works by tagging contigs based on their taxonomic classification (via Kraken 2 and TaxonKit) and examining their linked, discrete membership components within the assembly graph. This allows clean, species-specific separation of contaminated assemblies. Would only recommend for "single" whole genome sequences. I have NOT tested this on metagenomic assemblies, but it may work!


## Installation
1. Install the software requirements
```bash
# Clone repository
git clone https://github.com/bananabenana/separate_contaminated_assembly_graphs
cd separate_contaminated_assembly_graphs

# Move to directory
cd separate_contaminated_assembly_graphs

# Use mamba (or optionally conda to install the required packages)
mamba env create -f environment.yml

# Activate environment
mamba activate separate_contaminated_assembly_graphs_env

# Test installation
python separate_contaminated_assembly_graphs.py -h
```
You need a few databases to use this tool.
2. Select a Kraken 2 database from here: https://github.com/BenLangmead/aws-indexes/blob/master/docs/k2.md
```bash
# For example, the 2025 PlusPFP Kraken 2 database:
wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20251015.tar.gz
tar -zxvf k2_pluspfp_20251015.tar.gz
```

3. Download taxdump from NCBI
```bash
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz
```

You're now ready to start de-contaminating your assemblies!


## Quick usage
Separate suspected contaminated assemblies at the `family` level - useful for distantly-related bacteria using 12 threads
```bash
python separate_contaminated_assembly_graphs.py \
  --gfa_directory assemblies/graphs \
  --fasta_dir assemblies/fastas \
  --output_dir separated_assemblies \
  --taxa_level "family" \
  --kraken2_db k2_pluspfp_20251015 \
  --taxonkit_db taxdump \
  --threads 12
```

## Outputs



## Requirements


## Citation
To come


## Authors
- Ben Vezina
