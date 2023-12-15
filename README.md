# Zuker-Implementation
Implementation of Zuker's algorithm in Python.
# How to run:
1. Clone or download the repository.
2. Add your fasta files.
3. `python3 zuker.py -f FASTA -d DBN [-l LIMIT]`
   
    a. FASTA: file name of fasta file
    
    b. DBN: file name of dbn file to output
    
    c. LIMIT (optional): limit the length of the sequence (The algorithm is not recommended for large inputs. Inputs over 1000 will take a while. Ex: a sequence of length >1400 took roughly 17 minutes when tested on a personal laptop).
5. Can run `python3 zuker.py --help` to see usage.
