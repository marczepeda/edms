# Endogenous Deep Mutational Scans (EDMS)
## Command Line Interface
```shell
edms -h # or edms <TAB>
```
## Package Organization
- gen: input/output, data wrangling, generating plots, and statistics.
    ```shell
    import edms.gen.cli as cli
    import edms.gen.image as im
    import edms.gen.io as io
    import edms.gen.plot as p
    import edms.gen.stat as st
    import edms.gen.tidy as t
    ```
- bio: molecular biology & tissue culture workflows.
    ```shell
    import edms.bio.clone as cl
    import edms.bio.fastq as fq
    import edms.bio.genbank as gb
    import edms.bio.ngs as ngs
    import edms.bio.pe as pe
    import edms.bio.pegLIT as pegLIT
    import edms.bio.primedesign as primedesign
    import edms.bio.qPCR as qPCR
    import edms.bio.sanger as sanger
    import edms.bio.signature as signature
    import edms.bio.transfect as tf
    ```
- dat: interacting with databases.
    ```shell
    import edms.dat.cosmic as co
    import edms.dat.cvar as cv
    import edms.dat.ncbi as ncbi
    ```

## PyPI Instructions
### Install
1. Install edms from PyPI
    ```shell
    pip install edms
    ```
2. Optional: fastq.py UMI methods need umi_tools, cutadapt, samtools, bowtie2, and fgbio in a seperate environment
    ```shell
    conda create -n umi_tools umi_tools cutadapt samtools bowtie2 fgbio
    ```

### Update
1. Update edms from PyPI
    ```shell
    pip install --upgrade edms
    ```

## GitHub Instructions
### Install
1. Download Anaconda:
    - Mac: https://docs.anaconda.com/anaconda/install/mac-os/
    - Windows: https://docs.anaconda.com/anaconda/install/windows/
    - Linux: https://docs.anaconda.com/anaconda/install/linux/
2. Download Git: https://github.com/git-guides/install-git
3. Clone edms from github:
    ```shell
    cd ~
    mkdir git
    cd git
    git clone https://github.com/marczepeda/edms.git
    cd edms 
    ```
4. Make the environment and install edms:
    ```shell
    conda env create -f edms.yml # When conda asks you to proceed, type "y"
    conda activate edms
    pip install -e . # Include the "."
    bash autocomplete.sh # Optional: follow CLI instructions
    conda deactivate
    ```
5. Optional: fastq.py UMI methods need umi_tools, cutadapt, samtools, bowtie2, and fgbio in a seperate environment
    ```shell
    conda create -n umi_tools umi_tools cutadapt samtools bowtie2 fgbio
    ```

### Update
1. Enter the environment and uninstall edms:
    ```shell
    cd ~/git/edms
    conda activate edms
    pip uninstall -y edms
    rm -rf build/ dist/ *.egg-info

    ```
2. Pull latest version from github and install edms:
    ```shell
    git pull origin main
    pip install -e . # Include the "."
    conda deactivate
    ```

## PE Strategies
| Strategy | Description | Reference |
|----------|-------------|---------- |
| PE1 | Cas9(H840A) - M-MLV RT<br>+ pegRNA | [Search-and-replace genome editing without double-strand breaks or donor DNA](https://www.nature.com/articles/s41586-019-1711-4) |
| PE2 | Cas9(H840A) – M-MLV RT(D200N/L603W/T330P/T306K/W313F)<br>+ pegRNA | [Search-and-replace genome editing without double-strand breaks or donor DNA](https://www.nature.com/articles/s41586-019-1711-4) |
| PE3 | Cas9(H840A) – M-MLV RT(D200N/L603W/T330P/T306K/W313F)<br>+ ngRNA (targets non-edited strand) | [Search-and-replace genome editing without double-strand breaks or donor DNA](https://www.nature.com/articles/s41586-019-1711-4) |
| PE4 | Cas9(H840A) – M-MLV RT(D200N/L603W/T330P/T306K/W313F)<br>+ MLH1dn (MMR evasion) | [Enhanced prime editing systems by manipulating cellular determinants of editing outcomes](https://www.cell.com/cell/fulltext/S0092-8674(21)01065-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421010655%3Fshowall%3Dtrue) |
| PE5 | Cas9(H840A) – M-MLV RT(D200N/L603W/T330P/T306K/W313F)<br>+ MLH1dn (MMR evasion)<br>+ ngRNA (targets non-edited strand) | [Enhanced prime editing systems by manipulating cellular determinants of editing outcomes](https://www.cell.com/cell/fulltext/S0092-8674(21)01065-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421010655%3Fshowall%3Dtrue) |
| PE6a-d | Cas9(H840A) – ...<br>PEa: ... - evo-Ec48 RT<br>PEb: ... - evo-Tf1 RT<br>PEc: ... - Tf1 RT variant<br>PEd: ... - M-MLV RT variant | [Phage-assisted evolution and protein engineering yield compact, efficient prime editors](https://www.cell.com/cell/fulltext/S0092-8674(23)00854-1?uuid=uuid%3Acdb9bfe9-fd83-4a51-8a65-51f2e8e5cfe2) |
| PE6e-f | Cas9(H840A) variants – ...<br>M-MLV RT(ΔRNAseH) | [Phage-assisted evolution and protein engineering yield compact, efficient prime editors](https://www.cell.com/cell/fulltext/S0092-8674(23)00854-1?uuid=uuid%3Acdb9bfe9-fd83-4a51-8a65-51f2e8e5cfe2) |
| PE7 | Cas9(H840A) – M-MLV RT(D200N/L603W/T330P/T306K/W313F) - La (RNA binding protein that stabilizes pegRNA)<br>+/- ngRNA (targets non-edited strand) | [Improving prime editing with an endogenous small RNA-binding protein](https://www.nature.com/articles/s41586-024-07259-6) |
| PEmax | Mammalian codon-optimized PE | [Enhanced prime editing systems by manipulating cellular determinants of editing outcomes](https://www.cell.com/cell/fulltext/S0092-8674(21)01065-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421010655%3Fshowall%3Dtrue) |
| pegRNA | spacer - scaffold - RTT - PBS (makes the edit) | [Search-and-replace genome editing without double-strand breaks or donor DNA](https://www.nature.com/articles/s41586-019-1711-4) |
| epegRNA | spacer - scaffold - RTT - PBS - linker - tevoPreQ (makes the edit; more stable pegRNA) | [Engineered pegRNAs improve prime editing efficiency](https://www.nature.com/articles/s41587-021-01039-7) |
| ngRNA | spacer - scaffold (targets non-edited strand) | [Search-and-replace genome editing without double-strand breaks or donor DNA](https://www.nature.com/articles/s41586-019-1711-4) |
| MLH1dn | Dominant negative MLH1 (MMR evasion) | [Enhanced prime editing systems by manipulating cellular determinants of editing outcomes](https://www.cell.com/cell/fulltext/S0092-8674(21)01065-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421010655%3Fshowall%3Dtrue) |
| silent mutations | Larger prime edits are more efficient through bypassing MMR | [Enhanced prime editing systems by manipulating cellular determinants of editing outcomes](https://www.cell.com/cell/fulltext/S0092-8674(21)01065-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421010655%3Fshowall%3Dtrue) |
| La | Small RNA binding protein that stabilizes pegRNA | [Improving prime editing with an endogenous small RNA-binding protein](https://www.nature.com/articles/s41586-024-07259-6) |
| PE-eVLP | Engineered Virus-Like Particle for Prime Editors | [Engineered virus-like particles for transient delivery of prime editor ribonucleoprotein complexes in vivo](https://www.nature.com/articles/s41587-023-02078-y) |
| dNTPs | HSCs have low dNTP levels, limiting reverse transcription | [Enhancing prime editing in hematopoietic stem and progenitor cells by modulating nucleotide metabolism](https://www.nature.com/articles/s41587-024-02266-4) |
| Vpx | HSCs express SAMHD1 (triphosphohydrolase), which depletes dNTPs. Accessory lentiviral protein Vpx, encoded by HIV-2 and simian immunodeficiency viruses (SIVs), associates with the CRL4-DCAF1 E3 ubiquitin ligase to target SAMHD1 for proteasomal degradation. | [Enhancing prime editing in hematopoietic stem and progenitor cells by modulating nucleotide metabolism](https://www.nature.com/articles/s41587-024-02266-4) |
| MLH-SB | Small protein binder that disrupts MLH1 & PMS2 binding (MMR evasion) | [AI-generated small binder improves prime editing (Preprint)](https://www.biorxiv.org/content/10.1101/2024.09.11.612443v1.full) |