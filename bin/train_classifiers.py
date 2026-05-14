# creating the taxonomic classifiers

## DADA2
import os

fasta_file    = "dna-sequences.fasta"
taxonomy_file = "taxonomy.tsv"
output_file   = "dna-sequences-with-taxonomy.fasta"

# load taxonomy into dict
taxonomy = {}
with open(taxonomy_file) as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            taxonomy[parts[0]] = parts[1]

# write combined fasta
with open(fasta_file) as f, open(output_file, "w") as out:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            seq_id = line[1:]
            taxon  = taxonomy.get(seq_id, "unclassified")
            out.write(f">{seq_id} {taxon}\n")
        else:
            out.write(line + "\n")

print(f"Done! Written to {output_file}")
            
## IDTAXA
import gzip
import shutil
import os

fasta_file    = "dna-sequences.fasta"
taxonomy_file = "taxonomy.tsv"
output_file   = "idtaxa-dna-sequences-with-taxonomy.fasta"

def format_taxonomy(taxon_string):
    parts = taxon_string.strip().split(";")
    cleaned = ["Root"]  # add Root at the beginning
    for part in parts:
        part = part.strip()
        if "__" in part:
            part = part.split("__", 1)[1]
        part = part.replace("[", "").replace("]", "").replace("_", " ")
        cleaned.append(part)
    return ";".join(cleaned) + ";"

# load taxonomy into dict
taxonomy = {}
with open(taxonomy_file) as f:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            taxonomy[parts[0]] = format_taxonomy(parts[1])

# write combined fasta
with open(fasta_file) as f, open(output_file, "w") as out:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            seq_id = line[1:]
            taxon  = taxonomy.get(seq_id, "Root;unclassified;")
            out.write(f">{seq_id} {taxon}\n")
        else:
            out.write(line + "\n")

# compress
with open(output_file, 'rb') as f_in:
    with gzip.open(output_file + '.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

os.remove(output_file)
print(f"Done! Written to {output_file}.gz")
