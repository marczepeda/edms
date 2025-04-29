# Endogenous Deep Mutational Scans (EDMS)
## Package Organization
- gen: input/output, data wrangling, generating plots, and statistics.
    ```shell
    import edms.gen.io as io
    import edms.gen.tidy as t
    import edms.gen.plot as p
    import edms.gen.stat as st
    import edms.gen.image as im
    import edms.gen.cli as cli
    ```
- bio: molecular biology & tissue culture workflows.
    ```shell
    import edms.bio.transfect as tf
    import edms.bio.ngs as ngs
    import edms.bio.fastq as fq
    import edms.bio.clone as cl
    import edms.bio.pe as pe
    import edms.bio.pegLIT as pegLIT
    import edms.bio.qPCR as qPCR
    import edms.bio.genbank as gb
    ```
- dat: interacting with databases.
    ```shell
    import edms.dat.cosmic as co
    import edms.dat.cvar as cv
    import edms.dat.ncbi as ncbi
    ```

## Examples
- See https://github.com/marczepeda/edms_ex

## Instructions
### Install
1. Download Anaconda.
    - Mac: https://docs.anaconda.com/anaconda/install/mac-os/
    - Windows: https://docs.anaconda.com/anaconda/install/windows/
    - Linux: https://docs.anaconda.com/anaconda/install/linux/
2. Download Git: https://github.com/git-guides/install-git
    - Check using command line terminal: git version
3. Download Docker. https://www.docker.com/
4. Make environment: write the following in a command line terminal...
    ```shell
    cd ~
    conda create --name edms python=3.12.2
    # When conda asks you to proceed, type "y" 
    
    conda activate edms
    mkdir git
    cd git
    ```
5. Download dependencies: write the following in a command line terminal...
    ```shell
    conda install pip
    conda install conda-forge::biopython
    
    pip install -U scikit-learn # Also, installs numpy, pandas, matplotlib, seaborn, scipy.
    
    conda install -c conda-forge statsmodels
    conda install anaconda::requests
    conda install bioconda::viennarna
    conda install conda-forge::python-levenshtein
    conda install conda-forge::adjusttext
    conda install -c conda-forge git
    pip install dna-features-viewer
    ```
6. Install edms: write the following in a command line terminal...
    ```shell
    git clone https://github.com/marczepeda/edms.git
    cd edms
    pip install -e .
    # Include the "."
    
    docker pull pinellolab/primedesign
    # Docker desktop app needs to be open
    
    conda deactivate
    ```
### Update
1. Enter environment & delete edms: write the following in a command line terminal...
    ```shell
    cd ~
    cd git
    conda activate edms
    pip uninstall edms
    # Enter 'Y' when prompted
    
    rm -r edms
    # Enter 'Y' three times to completely remove the folder
    ```
2. Install edms: write the following in a command line terminal...
    ```shell
    git clone https://github.com/marczepeda/edms.git
    cd edms
    pip install -e .
    # Include the "."

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