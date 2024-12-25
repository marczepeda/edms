''' 
Module: ncbi.py
Author: Marc Zepeda
Created: 2024-12-25
Description: National Center of Biotechnology Information

Usage:
[]
- 
'''
# Import packages
from Bio import Entrez

'''
Entrez.email = "your_email@example.com"
handle = Entrez.esearch(db="protein", term="hemoglobin", retmax=10)
record = Entrez.read(handle)
print(record["IdList"])
'''