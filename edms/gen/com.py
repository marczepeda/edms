''' 
Module: com.py
Author: Marc Zepeda
Created: 2025-4-01
Description: Command Line Interaction

Usage:
[Common commands]
- expand_subs(): delete subdirectories and move their files to the parent directory

'''

# Import packages
import subprocess

# Common commands
def expand_subs(pt: str):
    ''' expand_subs(): delete subdirectories and move their files to the parent directory
        
        Parameters:
        pt (str): path to parent directory

        Dependencies: subprocess
    '''
    command = 'find . -not -type d -print0 | xargs -0J % mv -f % . ; find . -type d -depth -print0 | xargs -0 rm -rf'
    print(f"Terminal\ncd {pt}\n{command}")
    result = subprocess.run(f"{command}", shell=True, cwd=pt, capture_output=True, text=True)
    print(f"Output:\n{result.stdout}")
    if result.stderr: print(f"Errors:\n{result.stderr}")