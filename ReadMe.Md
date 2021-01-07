### Prereqs
- Python3.8
- pip3.8
- virtualenv
- a sample fasta file (e.g. dna_files/sample.fasta)

### Setup
1. `mkdir venv`
1. `virtualenv venv -p python3.8`
1. `pip3.8 install -r requirements.txt`
1. `make run`


### Result Format
- Results are in `dna_files/results`
- Each result file follows this naming convention:
  - `{input_name}.d{dist}_c{count_of_top_unique_seq}_g{num_groups}_t{total_input_seq}.fasta`
  - e.g. `sample.d1_c740_g40_t805.fasta`

### Notes (Updated: Dec 20, 2019)
- Each string contains 4 letters, and letters can only be A, T, C, G.
- Each string is ~200 letters.
- Each file contains ~10k strings.
- Similarity is dependent on the position of the letter, and not just based of # of similar letters.
- If 190/200 of a string of 200 characters are the same, then we consider them as the same.

### Requirements
- Read all strings from a fasta file, and find unique strings (don't need to compare strings across different files).
- Find similar strings, and group them together (e.g. AATT and AAAT, etc).
- Use % for similarity.

### Sample Output
```
> sample-23 seq appeared 740 times
aaccgg
> sample-25 seq appeared 30 times
aacggg
> sample-80 seq appeared 3 times
ccaagg
...
```
