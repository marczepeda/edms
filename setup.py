from setuptools import find_packages, setup
import os
work_dir = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(work_dir,'README.md'),'r') as f:
    long_description = f.read()

setup(
    name='edms',
    version='0.0.1',
    packages=find_packages(),
    description='Endogenous Deep Mutational Scans',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/marczepeda/edms',
    author='Marcanthony Zepeda',
    author_email='mzepeda@g.harvard.edu',
    license='MIT',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=['numpy',
                      'matplotlib',
                      'seaborn',
                      'pandas',
                      'biopython',
                      'requests',
                      'scipy',
                      'statsmodels',
                      'scikit-learn',
                      'Levenshtein',
                      'ViennaRNA',
                      'adjustText',
                      'pillow'],
    python_requires='>=3.11'
)