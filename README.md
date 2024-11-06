# Endogenous Deep Mutational Scans (EDMS)
## Package Organization
- gen: input/output, data wrangling, generating plots, and statistics.
    - import edms.gen.io as io
    - import edms.gen.tidy as t
    - import edms.gen.plot as p
    - import edms.gen.stat as st
- bio: molecular biology workflows including cloning and sequencing.
    - import edms.bio.fastq as f
    - import edms.bio.clone as cl
    - import edms.bio.pe as pe
    - import edms.bio.pegLIT as pegLIT
- dat: interacting with databases.
    - import edms.dat.cosmic as co
    - import edms.dat.cvar as cv

## Examples
- Design Library Workflow: FOXA1 a3_wing2 Domain
    - See /examples/design_library/src/foxa1_a3_wing2.ipynb
- Analyze Fastq Workflow: IKZF1 ZF2,3 Domain
    - See /examples/analyze_fastq/src/ikzf1_zf2_3.ipynb

## Instructions
### Clone (Option 1)
1. Download Anaconda.
    - Mac: https://docs.anaconda.com/anaconda/install/mac-os/
    - Windows: https://docs.anaconda.com/anaconda/install/windows/
    - Linux: https://docs.anaconda.com/anaconda/install/linux/
2. Download Git: https://github.com/git-guides/install-git
    - Check using command line terminal: git version
3. Download Docker. https://www.docker.com/
4. Write the following in a command line terminal:
    - cd ~
    - conda create --name emds_env python=3.11.5
        - When conda asks you to proceed, type "y" 
    - conda activate edms_env
    - mkdir git_edms
    - cd git_edms
    - git clone https://github.com/marczepeda/edms.git
    - pip install -e .
    - docker pull pinellolab/primedesign
        - Docker desktop app needs to be open

### Install (Option 2)
1. Download Anaconda.
    - Mac: https://docs.anaconda.com/anaconda/install/mac-os/
    - Windows: https://docs.anaconda.com/anaconda/install/windows/
    - Linux: https://docs.anaconda.com/anaconda/install/linux/
2. Download Git: https://github.com/git-guides/install-git
    - Check using command line terminal: git version
3. Download Docker: https://www.docker.com/
4. Write the following in a command line terminal:
    - conda create --name edms_env python=3.11.5
        - When conda asks you to proceed, type "y" 
    - conda activate edms_env
    - pip install git+https://github.com/marczepeda/edms.git
        - Sign into github when asked.
        - Wait at least a minute for the authentication to process.
    - conda list
        - check for edms package
    - docker pull pinellolab/primedesign
        - Docker desktop app needs to be open

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