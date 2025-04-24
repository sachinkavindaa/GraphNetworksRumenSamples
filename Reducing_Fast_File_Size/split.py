from pathlib import Path

def split_fasta(input_file, output_prefix, chunk_size=1000000):
    input_path = Path(input_file)
    output_dir = Path(output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    with input_path.open('r') as infile:
        chunk_number = 1
        current_chunk = []
        current_size = 0

        for line in infile:
            if line.startswith(">"):  # Start of a new sequence
                if current_chunk and current_size >= chunk_size:
                    output_file = output_dir / f"{Path(output_prefix).stem}_chunk_{chunk_number}.faa"
                    with output_file.open('w') as chunk_file:
                        chunk_file.writelines(current_chunk)
                    print(f"✅ Wrote chunk {chunk_number}: {len(current_chunk)} lines")
                    chunk_number += 1
                    current_chunk = []
                    current_size = 0
            current_chunk.append(line)
            current_size += len(line)

        # Write the last chunk
        if current_chunk:
            output_file = output_dir / f"{Path(output_prefix).stem}_chunk_{chunk_number}.faa"
            with output_file.open('w') as chunk_file:
                chunk_file.writelines(current_chunk)
            print(f"✅ Wrote final chunk {chunk_number}: {len(current_chunk)} lines")

    print(f"\n✅ FASTA file split into {chunk_number} chunks in {output_dir}")

# Call the function
split_fasta(
    "/work/samodha/sachin/GraphN/dedup_filtered.faa",
    "/work/samodha/sachin/GraphN/split_fasta/dedup_filtered_split"
)
