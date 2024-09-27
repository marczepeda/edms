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
- dat: interacting with databases
    - import edms.dat.cosmic as co
    - import edms.dat.cvar as cv

## Examples
- Design Library Workflow: FOXA1 a3_wing2 Domain
    - See /examples/design_library/src/foxa1_a3_wing2.ipynb
- Analyze Fastq Workflow: IKZF1 ZF2,3 Domain
    - See /examples/analyze_fastq/src/ikzf1_zf2_3.ipynb

## Instructions
### Install (Option 1)
1. Open Terminal.
2. Write the following commands:
    - conda create --name emds_env
        - When conda asks you to proceed, type "y" 
    - conda activate edms_env
    - pip install git+https://github.com/marczepeda/edms.git

### Clone (Option 2)
1. Open Terminal.
2. Write the following commands:
    - cd
    - conda create --name emds_env
        - When conda asks you to proceed, type "y" 
    - conda activate edms_env
    - mkdir edms_git
    - cd edms_git
    - git clone https://github.com/marczepeda/edms.git
    - pip -e .

## PE Strategies
| Strategy | Description | Reference |
|----------|-------------|-----------|
| PE1 | Cas9(H840A) - M-MLV RT + pegRNA | [Search-and-replace genome editing without double-strand breaks or donor DNA](https://www.nature.com/articles/s41586-019-1711-4) |
| PE2 | Cas9(H840A) – M-MLV RT(D200N/L603W/T330P/T306K/W313F) + pegRNA | [Search-and-replace genome editing without double-strand breaks or donor DNA](https://www.nature.com/articles/s41586-019-1711-4) |
| PE3 | PE2 + additional ngRNA (targets non-edited strand) | [Search-and-replace genome editing without double-strand breaks or donor DNA](https://www.nature.com/articles/s41586-019-1711-4) |
| PE4 | PE2 + MLH1dn (MMR evasion) | [Enhanced prime editing systems by manipulating cellular determinants of editing outcomes](https://www.cell.com/cell/fulltext/S0092-8674(21)01065-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421010655%3Fshowall%3Dtrue) |
| PE5 | PE2 + MLH1dn (MMR evasion) + additional ngRNA (targets non-edited strand) | [Enhanced prime editing systems by manipulating cellular determinants of editing outcomes](https://www.cell.com/cell/fulltext/S0092-8674(21)01065-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421010655%3Fshowall%3Dtrue) |
| PE6 | PE2 - La (RNA binding protein that stabilizes pegRNA) | [Improving prime editing with an endogenous small RNA-binding protein](https://www.nature.com/articles/s41586-024-07259-6) |
| PE7 | PE2 - La (RNA binding protein that stabilizes pegRNA) + (targets non-edited strand) | [Improving prime editing with an endogenous small RNA-binding protein](https://www.nature.com/articles/s41586-024-07259-6) |
| PEmax | Mammalian codon-optimized PE | [Enhanced prime editing systems by manipulating cellular determinants of editing outcomes](https://www.cell.com/cell/fulltext/S0092-8674(21)01065-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421010655%3Fshowall%3Dtrue) |
| pegRNA | spacer - scaffold - RTT - PBS (makes the edit) | [Search-and-replace genome editing without double-strand breaks or donor DNA](https://www.nature.com/articles/s41586-019-1711-4) |
| epegRNA | spacer - scaffold - RTT - PBS - linker - tevoPreQ (makes the edit; more stable) | [Engineered pegRNAs improve prime editing efficiency](https://www.nature.com/articles/s41587-021-01039-7) |
| ngRNA | spacer - scaffold (targets non-edited strand) | [Search-and-replace genome editing without double-strand breaks or donor DNA](https://www.nature.com/articles/s41586-019-1711-4) |
| MLH1dn | Dominant negative MLH1 (MMR evasion) | [Enhanced prime editing systems by manipulating cellular determinants of editing outcomes](https://www.cell.com/cell/fulltext/S0092-8674(21)01065-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421010655%3Fshowall%3Dtrue) |
| silent mutations | Larger prime edits are more efficient through bypassing MMR | [Enhanced prime editing systems by manipulating cellular determinants of editing outcomes](https://www.cell.com/cell/fulltext/S0092-8674(21)01065-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421010655%3Fshowall%3Dtrue) |
| La | Small RNA binding protein that stabilizes pegRNA | [Improving prime editing with an endogenous small RNA-binding protein](https://www.nature.com/articles/s41586-024-07259-6) |
| PE-eVLP | Engineered Virus-Like Particle for Prime Editors | [Engineered virus-like particles for transient delivery of prime editor ribonucleoprotein complexes in vivo](https://www.nature.com/articles/s41587-023-02078-y) |
| dNTP supplementation | HSCs have low dNTP levels, limiting reverse transcription | [Enhancing prime editing in hematopoietic stem and progenitor cells by modulating nucleotide metabolism](https://www.nature.com/articles/s41587-024-02266-4) |
| VPX | HSCs express SAMHD1 (triphosphohydrolase), which depletes dNTPs. Accessory lentiviral protein Vpx, encoded by HIV-2 and simian immunodeficiency viruses (SIVs), associates with the CRL4-DCAF1 E3 ubiquitin ligase to target SAMHD1 for proteasomal degradation | [Enhancing prime editing in hematopoietic stem and progenitor cells by modulating nucleotide metabolism](https://www.nature.com/articles/s41587-024-02266-4) |
| MLH-SB | Small protein binder that disrupts MLH1 & PMS2 binding (MMR evasion) | [AI-generated small binder improves prime editing (Preprint)](https://www.biorxiv.org/content/10.1101/2024.09.11.612443v1.full) |