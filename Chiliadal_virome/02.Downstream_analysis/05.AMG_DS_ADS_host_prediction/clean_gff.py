import sys

if len(sys.argv) != 3:
    print("Usage: python clean_gff.py <input.gff> <output.gff>")
    sys.exit(1)

input_gff = sys.argv[1]
output_gff = sys.argv[2]

id_counts = {}

with open(input_gff, "r") as infile, open(output_gff, "w") as outfile:
    for line in infile:
        line = line.strip()
        if not line or line.startswith("##gff-version"):
            continue  # skip version headers and blank lines

        fields = line.split("\t")
        if len(fields) < 9:
            continue  # skip malformed lines

        seqid = fields[0]
        attributes = fields[8].rstrip(";")  # remove trailing semicolon if present

        # Increment ID counter for this seqid
        id_counts[seqid] = id_counts.get(seqid, 0) + 1
        new_id = f"ID={seqid}_{id_counts[seqid]}"
        
        # --- process the attributes ---
        if attributes.startswith("ID="):
            # Split off the existing ID
            parts = attributes.split(";", 1)
            old_id_value = parts[0].replace("ID=", "", 1)  # remove the "ID=" prefix
            remaining = parts[1] if len(parts) > 1 else ""

            # Construct new attribute order
            if remaining:
                new_attributes = f"{new_id};{old_id_value};{remaining}"
            else:
                new_attributes = f"{new_id};{old_id_value}"
        else:
            # No existing ID, just prepend the new one
            if attributes:
                new_attributes = f"{new_id};{attributes}"
            else:
                new_attributes = new_id

        fields[8] = new_attributes + ";"
        outfile.write("\t".join(fields) + "\n")
