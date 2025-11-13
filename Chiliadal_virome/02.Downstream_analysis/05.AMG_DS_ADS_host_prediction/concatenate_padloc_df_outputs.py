#######################################################
# Script was created by dr. James Docherty
#######################################################

import pandas as pd
import sys

# Check for correct usage
if len(sys.argv) != 4:
    print("Usage: python convert_padloc_to_DF.py <padloc_file> <defensefinder_file> <combined_output>")
    sys.exit(1)

padloc_file = sys.argv[1]
defensefinder_file = sys.argv[2]
combined_output = sys.argv[3]

# read both files
padloc_df = pd.read_csv(padloc_file, sep="\t", dtype=str)
defense_df = pd.read_csv(defensefinder_file, sep="\t", dtype=str)

# recover seqid in defensefinder from sys_id
def get_seqid(sys_beg):
    return sys_beg.rpartition("_")[0]

defense_df["seqid"] = defense_df["sys_beg"].apply(get_seqid)

# add missing 'activity' to padloc if needed
if 'activity' not in padloc_df.columns:
    padloc_df['activity'] = 'Defense'

# add source labels
padloc_df["source"] = "padloc"
defense_df["source"] = "defensefinder"

# define consistent columns
desired_columns = [
    "sys_id", "type", "subtype", "activity", "sys_beg", "sys_end",
    "protein_in_syst", "genes_count", "name_of_profiles_in_sys", "seqid", "source"
]

# ensure columns match
for col in desired_columns:
    if col not in padloc_df.columns:
        padloc_df[col] = ''
    if col not in defense_df.columns:
        defense_df[col] = ''

padloc_df = padloc_df[desired_columns]
defense_df = defense_df[desired_columns]

# combine them
combined_df = pd.concat([defense_df, padloc_df], ignore_index=True)

# remove fully empty rows
combined_df.dropna(how="all", inplace=True)

# save to output
combined_df.to_csv(combined_output, sep="\t", index=False)

print(f"✅ Combined file written to {combined_output}")

