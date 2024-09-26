# Endogenous Deep Mutational Scans (EDMS)
## Package Organization
- gen: input/output, data wrangling, generating plots, and statistics.\
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