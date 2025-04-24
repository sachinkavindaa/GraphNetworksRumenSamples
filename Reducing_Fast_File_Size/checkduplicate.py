import hashlib
from collections import defaultdict

fasta_file = "/work/samodha/sachin/GraphN/filtered.faa"
hash_counts = defaultdict(int)

def hash_seq(seq):
    return hashlib.md5(seq.encode()).hexdigest()

with open(fasta_file, "r") as f:
    current_seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                seq_hash = hash_seq(current_seq)
                hash_counts[seq_hash] += 1
                current_seq = ""
        else:
            current_seq += line.upper()
    if current_seq:
        seq_hash = hash_seq(current_seq)
        hash_counts[seq_hash] += 1

duplicates = {h: c for h, c in hash_counts.items() if c > 1}

print(f"Total unique sequences (by hash): {len(hash_counts)}")
print(f"Number of duplicated sequences (by hash): {len(duplicates)}")

if duplicates:
    print("\nExamples of duplicated hashes:")
    for i, (seq_hash, count) in enumerate(duplicates.items()):
        print(f"Duplicate #{i+1}: {count} copies — Hash: {seq_hash}")
        if i >= 4:
            break
