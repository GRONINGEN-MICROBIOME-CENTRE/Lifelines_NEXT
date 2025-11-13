#######################################################
# Script was created by dr. James Docherty
#######################################################

import pandas as pd
import sys

# Check for correct usage
if len(sys.argv) != 3:
    print("Usage: python convert_padloc_to_DF.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    # read combined tab-separated
    df = pd.read_csv(input_file, sep="\t", dtype=str)

    # filter out PDC systems and CRISPR_array
    df = df[~df['system'].str.startswith('PDC', na=False)]
    df = df[df['system'] != 'CRISPR_array']

    # group by genome + system.number
    grouped = df.groupby(['seqid', 'system.number'])

    # DefenseFinder-style columns
    result = pd.DataFrame(columns=[
        'sys_id', 'type', 'subtype', 'sys_beg', 'sys_end',
        'protein_in_syst', 'genes_count', 'name_of_profiles_in_sys', 'seqid'
    ])

    for (seqid, sys_num), group in grouped:
        new_row = pd.DataFrame({
            'sys_id': [f"{seqid}_{sys_num}"],
            'type': [''],  # blank for now
            'subtype': [group['system'].iloc[0]],
            'sys_beg': [group['target.name'].iloc[0]],
            'sys_end': [group['target.name'].iloc[-1]],
            'protein_in_syst': [','.join(group['target.name'])],
            'genes_count': [len(group)],
            'name_of_profiles_in_sys': [','.join(group['hmm.name'].dropna().unique())],
            'seqid': [seqid]
        })
        result = pd.concat([result, new_row], ignore_index=True)

    # sort and reset
    result = result.sort_values(['seqid', 'sys_id']).reset_index(drop=True)

    # save
    result.to_csv(output_file, sep="\t", index=False)

    print(f"✅ Parsed DefenseFinder-style summary written to {output_file}")

except Exception as e:
    print(f"❌ Error: {e}")

