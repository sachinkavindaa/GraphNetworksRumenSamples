import hashlib

input_fasta = "/work/samodha/sachin/GraphN/filtered.faa"
dedup_fasta = "/work/samodha/sachin/GraphN/dedup_filtered.faa"
seen_hashes = set()

def hash_seq(seq):
    return hashlib.md5(seq.encode()).hexdigest()

with open(input_fasta, "r") as infile, open(dedup_fasta, "w") as outfile:
    current_header = ""
    current_seq = ""

    for line in infile:
        line = line.strip()
        if line.startswith(">"):
            if current_header and current_seq:
                h = hash_seq(current_seq.upper())
                if h not in seen_hashes:
                    seen_hashes.add(h)
                    outfile.write(current_header + "\n")
                    for i in range(0, len(current_seq), 60):
                        outfile.write(current_seq[i:i+60] + "\n")
                current_seq = ""
            current_header = line
        else:
            current_seq += line

    # Handle last sequence
    if current_header and current_seq:
        h = hash_seq(current_seq.upper())
        if h not in seen_hashes:
            seen_hashes.add(h)
            outfile.write(current_header + "\n")
            for i in range(0, len(current_seq), 60):
                outfile.write(current_seq[i:i+60] + "\n")

print(f"Finished writing deduplicated FASTA to: {dedup_fasta}")
print(f"Total unique sequences written: {len(seen_hashes)}")
