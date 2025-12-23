#!/bin/bash
#SBATCH --job-name=PDB_Uniprot_annotations
#SBATCH --output=pdb_uniprot_annotations.out
#SBATCH --mem=20gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --open-mode=truncate
#SBATCH --partition=regular

IDs=$1       # file with PDB IDs
output=$2    # output file

# Header
echo -e "PDB_ID\tPDB_title\tUniProt_ID\tProtein_name\tFunction\tGO_CC\tGO_MF\tGO_BP\tKeywords\tDomains" > "$output"

while read -r PDBID; do
    [ -z "$PDBID" ] && continue
    echo "Processing $PDBID..." >&2

    # Get PDB title
    PDB_TITLE=$(curl -s "https://data.rcsb.org/rest/v1/core/entry/$PDBID" | jq -r '.struct.title // "-"')

    # Submit mapping job: PDB -> UniProt
    JOB_ID=$(curl -sf -X POST "https://rest.uniprot.org/idmapping/run" \
        -d from=PDB -d to=UniProtKB -d ids="$PDBID" | jq -r '.jobId')

    if [ -z "$JOB_ID" ] || [ "$JOB_ID" == "null" ]; then
        echo "Failed to submit mapping for $PDBID" >&2
        continue
    fi

    # Wait for job to finish
    STATUS="RUNNING"
    COUNT=0
    while [ "$STATUS" != "FINISHED" ] && [ $COUNT -lt 60 ]; do
        sleep 2
        STATUS=$(curl -sf "https://rest.uniprot.org/idmapping/status/$JOB_ID" | jq -r '.jobStatus // "RUNNING"')
        COUNT=$((COUNT+1))
    done

    if [ "$STATUS" != "FINISHED" ]; then
        echo "Timeout waiting for mapping of $PDBID" >&2
        continue
    fi

    # Get mapped UniProt IDs
    UNI_IDS=$(curl -sf "https://rest.uniprot.org/idmapping/stream/$JOB_ID" | jq -r '.results[].to // empty')

    if [ -z "$UNI_IDS" ]; then
        echo "No UniProt ID found for $PDBID" >&2
        continue
    fi

    # Retrieve UniProt info for each UniProt ID
    for UNIPROT in $UNI_IDS; do
        echo "Fetching annotations for $UNIPROT..." >&2

        curl -s "https://rest.uniprot.org/uniprotkb/${UNIPROT}.json" | jq -r --arg pdb "$PDBID" --arg title "$PDB_TITLE" '
        def safe_join($arr): ($arr | join("; ")) as $res | if $res=="" then "-" else $res end;
        def safe_value($v): if $v == null or $v=="" then "-" else $v end;

        [
            $pdb,
            $title,
            safe_value(.primaryAccession),
            safe_value(.proteinDescription.recommendedName.fullName.value),

            # Function: first sentence, sanitized
            (
              (
                [.comments[]? | select(.commentType=="FUNCTION") | .texts[]?.value]
                | join(" ")
              )
              | gsub("[\t\n\r]+";" ")
              | capture("^(?<first>.*?\\.)")?
              | safe_value(.first)
            ),

            # GO CC
            safe_join([.uniProtKBCrossReferences[]? | select(.database=="GO" and (.properties[0].value | startswith("C"))) | "\(.id):\([.properties[]?.value] | join("; "))"]),

            # GO MF
            safe_join([.uniProtKBCrossReferences[]? | select(.database=="GO" and (.properties[0].value | startswith("F"))) | "\(.id):\([.properties[]?.value] | join("; "))"]),

            # GO BP
            safe_join([.uniProtKBCrossReferences[]? | select(.database=="GO" and (.properties[0].value | startswith("P"))) | "\(.id):\([.properties[]?.value] | join("; "))"]),

	    # Pathway DBs
     	    safe_value([.uniProtKBCrossReferences[]? | select(.database=="Reactome" or .database=="PathwayCommons" or .database=="SIGNOR" or .database=="SignaLink") | .database] | join(", ")),

	    # Pathway IDs
	    safe_value([.uniProtKBCrossReferences[]? | select(.database=="Reactome" or .database=="PathwayCommons" or .database=="SIGNOR" or .database=="SignaLink") | .id] | join(", ")),
 
	    # Pathway descriptions
	    safe_value([.uniProtKBCrossReferences[]? | select(.database=="Reactome" or .database=="PathwayCommons" or .database=="SIGNOR" or .database=="SignaLink") | (.properties[0].value // "-")] | join("; ")),

            # Keywords
            safe_join([.keywords[]? | .name]),

            # Domains
            safe_join([.uniProtKBCrossReferences[]? | select(.database=="Pfam" or .database=="InterPro" or .database=="SUPFAM" or .database=="Gene3D") | "\(.database):\(.id)"])
        ] | @tsv
        ' >> "$output"

    done
done < "$IDs"

echo "Annotation completed. Output written to $output" >&2
