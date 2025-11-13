import sys

# Check for correct usage
if len(sys.argv) != 3:
    print("Usage: python clean_faa.py <input.faa> <output.faa>")
    sys.exit(1)

input_faa = sys.argv[1]
output_faa = sys.argv[2]

id_counts = {}

with open(input_faa, "r") as infile, open(output_faa, "w") as outfile:
    current_seq_lines = []
    current_header = ""

    for line in infile:
        line = line.strip()
        if line.startswith(">"):
            # Write the previous entry, if any
            if current_header:
                outfile.write(f"{current_header}\n")
                outfile.write("\n".join(current_seq_lines) + "\n")

            # Strip '>' and remove everything from "_CDS_[" onward
            base_id = line[1:].split("_CDS_[")[0]
            id_counts[base_id] = id_counts.get(base_id, 0) + 1
            new_header = f">{base_id}_{id_counts[base_id]}"
            current_header = new_header
            current_seq_lines = []
        else:
            current_seq_lines.append(line)

    # Write the final entry
    if current_header:
        outfile.write(f"{current_header}\n")
        outfile.write("\n".join(current_seq_lines) + "\n")

