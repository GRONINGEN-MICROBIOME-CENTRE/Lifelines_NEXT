#!/bin/bash
#SBATCH --job-name=Uniprot_annotations
#SBATCH --output=Uniprot_annotations.out
#SBATCH --mem=20gb
#SBATCH --time=00-2:00 
#SBATCH --cpus-per-task=32
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

IDs=$1   # file with UniProt IDs
output=$2

# Add header
echo -e "UniProt_ID\tProtein_name\tFunction_short\tGO_CC\tGO_MF\tGO_BP\tPathway_DBs\tPathway_IDs\tPathway_descriptions\tKeywords\tDomains" > "$output"

while read -r ID; do
  [ -z "$ID" ] && continue
  echo "Fetching $ID..." >&2

  curl -s "https://rest.uniprot.org/uniprotkb/${ID}.json" | jq -r '
    # Function to replace empty/null with -
    def safe_value($val): if $val == "" or $val == null then "-" else $val end;

    [
      (.primaryAccession // "-"),
      (.proteinDescription.recommendedName.fullName.value // "-"),

      # Function: first sentence, fully safe
      (
        try (
          [.comments[]? | select(.commentType=="FUNCTION") | .texts[]?.value]
          | map(tostring)
          | join(" ")
          | gsub("[\t\n\r]+";" ")
          | capture("^(?<first>.*?\\.)")?
          | .first
        ) catch "-"
      ) // "-",

      # GO CC
      safe_value([.uniProtKBCrossReferences[]? | select(.database=="GO" and (.properties[0].value | startswith("C"))) | "\(.id):\([.properties[]?.value] | join("; "))"] | join(", ")),

      # GO MF
      safe_value([.uniProtKBCrossReferences[]? | select(.database=="GO" and (.properties[0].value | startswith("F"))) | "\(.id):\([.properties[]?.value] | join("; "))"] | join(", ")),

      # GO BP
      safe_value([.uniProtKBCrossReferences[]? | select(.database=="GO" and (.properties[0].value | startswith("P"))) | "\(.id):\([.properties[]?.value] | join("; "))"] | join(", ")),

      # Pathway DBs
      safe_value([.uniProtKBCrossReferences[]? | select(.database=="Reactome" or .database=="PathwayCommons" or .database=="SIGNOR" or .database=="SignaLink") | .database] | join(", ")),

      # Pathway IDs
      safe_value([.uniProtKBCrossReferences[]? | select(.database=="Reactome" or .database=="PathwayCommons" or .database=="SIGNOR" or .database=="SignaLink") | .id] | join(", ")),

      # Pathway descriptions
      safe_value([.uniProtKBCrossReferences[]? | select(.database=="Reactome" or .database=="PathwayCommons" or .database=="SIGNOR" or .database=="SignaLink") | (.properties[0].value // "-")] | join("; ")),

      # Keywords
      safe_value([.keywords[]? | .name] | join("; ")),

      # Domains
      safe_value([.uniProtKBCrossReferences[]? | select(.database=="InterPro" or .database=="Pfam" or .database=="Gene3D" or .database=="SUPFAM") | "\(.database):\(.id)"] | join("; "))
    ] | @tsv
  ' >> "$output"

done < "$IDs"

echo "Annotation completed. Output written to $output" >&2

