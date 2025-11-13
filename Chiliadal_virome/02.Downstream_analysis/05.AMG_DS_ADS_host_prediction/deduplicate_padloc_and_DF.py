#######################################################
# Script was created by dr. James Docherty
#######################################################

import pandas as pd
import sys

# Check for correct usage
if len(sys.argv) != 5:
    print("Usage: python convert_padloc_to_DF.py <input_file> <deduplicated_out> <successfully_deduplicated_out> <non_deduplicated_out>")
    sys.exit(1)

input_file = sys.argv[1]
deduplicated_out = sys.argv[2]
successfully_deduplicated_out = sys.argv[3]
non_deduplicated_out = sys.argv[4]

# Load the combined file
df = pd.read_csv(input_file, sep="\t", dtype=str)

# Helper function for tool_final logic
def get_tool_final(subtype, tools, original_tool):
    if subtype == "mza":
        return "padloc"
    elif subtype in {"Cas_Type_II-C", "Abi2"}:
        return original_tool
    elif len(tools) == 2:
        return "both"
    return tools.pop()

# Rules for merging conflicting subtypes
subtype_rules = {
    ("Lamassu_Family", "Lamassu-Mrr"): "Lamassu-Mrr",
    ("RM_IIG", "BREX_I"): "BREX_I",
    ("RM_II", "RM_HNH"): "RM_HNH",
    ("RM_I", "PrrC"): "PrrC",
    ("Cas_Type_II-C", "Cas_Cluster"): "Cas_Type_II-C",
    ("mza", "RM_II"): "mza",
    ("Septu", "Septu_I"): "Septu_I",
    ("Thoeris_I", "Thoeris_II"): "Thoeris_II",
    ("AbiD", "Abi2"): "Abi2",
    ("Hachiman_I", "Hachiman"): "Hachiman_I",
    ("Cas_Type_II-C", "Cas_Type_II-A"): "Cas_Type_II-C",
    ("kiwa", "Kiwa"): "Kiwa",
    ("RM_type_IIG", "RM_Type_IIG"): "RM_Type_IIG",
    ("DRT_class_II", "UG27"): "UG27",
    ("dXTPase", "Ogmios"): "Ogmios",
    ("RM_type_II", "RM_Type_II"): "RM_Type_II",
    ("retron_II-A", "Retron_II"): "Retron_II-A",
    ("PD-T7-1or5", "PD-T7-1"): "PD-T7-1",
    ("Lamassu_Family", "Lamassu-Hydrolase_Protease"): "Lamassu_Family",
    ("DRT_other", "Retron_VI"): "Retron_VI",
    ("Paris_fused", "PARIS_II_merge"): "PARIS_II_merge",
    ("rexb", "RexAB"): "RexAB",
    ("darTG", "DarTG"): "DarTG",
    ("septu_type_I", "Septu"): "Septu_type_I",
    ("RM_type_IV", "RM_Type_IV"): "RM_Type_IV",
    ("retron_I-C", "Retron_I_C"): "Retron_I_C",
    ("shedu", "Shedu"): "Shedu",
    ("GAO_19", "PD-T7-2"): "PD-T7-2",
    ("dsr2", "Dsr_II"): "Dsr_II",
    ("gabija", "Gabija"): "Gabija",
    ("hachiman_type_I", "Hachiman"): "Hachiman_type_I",
    ("pycsar_effector", "Pycsar"): "Pycsar",
    ("cbass_type_I", "CBASS_I"): "CBASS_I",
    ("wadjet_other", "Wadjet_III"): "Wadjet_III",
    ("RM_type_I", "RM_Type_I"): "RM_Type_I",
    ("Lamassu_Family", "Lamassu-Cap4_nuclease"): "Lamassu_Family",
    ("Paris", "PARIS_II"): "PARIS_II",
    ("retron_I-B", "Retron_I_B"): "Retron_I_B",
    ("gop_beta_cll", "Rst_gop_beta_cll"): "Rst_gop_beta_cll",
    ("Lamassu_Family", "Lamassu-Hypothetical"): "Lamassu_Family",
    ("retron_II-A", "Retron_II-A"): "Retron_II-A",
    ("septu_type_I", "Septu_type_I"): "Septu_type_I",
    ("DMS_other", "RM_Type_I"): "RM_Type_I",
    ("DMS_other", "RM_Type_III"): "RM_Type_III",
    ("RM_type_III", "RM_Type_III"): "RM_Type_III",
    ("DMS_other", "RM_Type_II"): "RM_Type_II",
    ("zorya_type_II", "Zorya_TypeII"): "Zorya_TypeII",
    ("dsr1", "Dsr_I"): "Dsr_I",
    ("cbass_type_IIs", "CBASS_II"): "CBASS_IIs",
    ("dXTPase", "dGTPase"): "dGTPase",
    ("acriia3", "acriia25"): "acriia",
    ("acriia25", "acriia3"): "acriia",
    ("hachiman_type_I", "Hachiman_type_I"): "Hachiman_type_I",
    ("cas_adaptation", "CAS_Class1-Subtype-I-E"): "CAS_Class1-Subtype-I-E",
    ("cas_type_other", "CAS_Class1-Subtype-I-A"): "CAS_Class1-Subtype-I-A",
    ("DMS_other", "hia5_hin1523_nma1821"): "hia5_hin1523_nma1821",
    ("DRT_class_I", "DRT8"): "DRT8",
    ("Helicase-DUF2290", "Rst_HelicaseDUF2290"): "Rst_HelicaseDUF2290",
    ("RM_type_IIG", "RM_Type_IIG_2"): "RM_Type_IIG_2",
    ("cas_type_IV-B", "CAS_Class1-Type-IV"): "CAS_Class1-Type-IV-B",
    ("cbass_type_IIs", "CBASS_IIs"): "CBASS_IIs",
    ("3HP", "Rst_3HP"): "Rst_3HP",
    ("argonaute_solo", "pAgo_S1A"): "pAgo_S1A",
    ("argonaute_solo", "pAgo_LongB"): "pAgo_LongB",
    ("ietAS", "Gao_Iet"): "Gao_Iet",
    ("cas_type_IV-B", "CAS_Class1-Subtype-IV-B"): "CAS_Class1-Subtype-IV-B",
    ("qatABCD", "Gao_Qat"): "Gao_Qat",
    ("VSPR", "RM_Type_II"): "RM_Type_II",
    ("cas_adaptation", "CAS_Cluster"): "CAS_Cluster",
    ("acriia25", "acriia"): "acriia25",
    ("acriia3", "acriia"): "acriia3",
    ("cas_type_IV-B", "CAS_Class1-Type-IV-B"): "CAS_Class1-Subtype-IV-B",
    ("RM_Type_II", "Shedu"): "Shedu",
    ("thoeris_type_I", "Thoeris_II"): "Thoeris_II",
    ("RM_type_II", "Shedu"): "Shedu",
    ("Old", "Old_exonuclease"): "Old",
    ("RM_type_IV", "HEC-06"): "HEC-06",
    ("RM_type_II", "RM_Type_IV"): "RM_type_II", 
    ("DRT_class_I", "AbiA_large"): "AbiA_large",
    ("AbiE", "SanaTA"): "SanaTA",
    ("Lamassu_Family", "RloC"): "RloC",
    ("cas_type_IV-B", "CAS_Class1-Subtype-IV-A"): "CAS_Class1-Subtype-IV-A",
    ("cas_type_IV-B", "CAS_Class1-Subtype-IV-D"): "CAS_Class1-Subtype-IV-D",
    ("Paris", "PARIS_I"): "PARIS_I",
    ("retron_IV", "Retron_IV"): "Retron_IV",
    ("retron_XII", "Retron_XII"): "Retron_XII",
    ("thoeris_type_I", "Thoeris_I"): "Thoeris_I",
    ("cas_type_I-C", "CAS_Class1-Subtype-I-C"): "CAS_Class1-Subtype-I-C",
    ("brex_type_IV", "BREX_IV"): "BREX_IV",
    ("retron_I-A", "Retron_I_A"): "Retron_I_A",
    ("Mokosh_TypeI", "Mokosh_Type_I_B"): "Mokosh_Type_I_B",
    ("druantia_other", "Druantia_I"): "Druantia_I",
    ("DMS_other", "dam"): "dam",
    ("DMS_other", "RM_Type_IV_1"): "RM_Type_IV_1",
    ("TIR-NLR", "Rst_TIR-NLR"): "Rst_TIR-NLR",
    ("AVAST_type_V", "Avs_V"): "AVAST_type_V",
    ("cas_type_other", "CAS_Class1-Type-IV-B"): "CAS_Class1-Type-IV-B",
    ("cas_type_II-C", "CAS_Cluster"): "CAS_Class1-Subtype-II-C",
    ("retron_III-A", "Retron_III"): "Retron_III-A",
    ("DMS_other", "PrrC"): "PrrC",
    ("cas_type_other", "CAS_Class1-Subtype-I-F"): "CAS_Class1-Subtype-I-F",
    ("cas_type_other", "CAS_Class1-Subtype-III-D"): "CAS_Class1-Subtype-III-D",
    ("Lamassu_Family", "Lamassu-Lipase"): "Lamassu_Family",
    ("mza", "Gao_Mza"): "mza",
    ("argonaute_type_II", "pAgo_LongB"): "pAgo_LongB",
    ("RM_type_I", "PrrC"): "PrrC",
    ("retron_III-A", "Retron_III-A"): "Retron_III-A",
    ("cas_type_II-C", "CAS_Class1-Subtype-II-C"): "CAS_Class1-Subtype-II-C",
    ("DRT_class_II", "UG17"): "UG17",
    ("DRT_class_III", "UG24"): "UG24",
    ("DRT_class_I", "AbiP2"): "AbiP2",
    ("DRT_class_I", "DRT_4"): "DRT_4",
    ("DRT_class_II", "DRT_2"): "DRT_2",
    ("dam", "RM_Type_II"): "RM_Type_II",
    ("RM_Type_II", "dam"): "RM_Type_II",
    ("Tin", "DS-18"): "Tin",
    ("GAO_19", "Gao_Her_SIR"): "PD-T7-2"
}

def resolve_subtype(subtype1, subtype2):
    if (subtype1, subtype2) in subtype_rules:
        return subtype_rules[(subtype1, subtype2)]
    elif (subtype2, subtype1) in subtype_rules:
        return subtype_rules[(subtype2, subtype1)]
    else:
        return None

# Containers
deduplicated = []
successfully_deduplicated = []
non_deduplicated = []

# group by genome (seqid)
grouped = df.groupby("seqid")

for genome, group in grouped:
    for index, row in group.iterrows():
        proteins = set(row["protein_in_syst"].split(","))
        overlap_found = False

        # try to match with existing deduplicated rows
        for dedup_row in deduplicated:
            if dedup_row["seqid"] == genome:
                existing_proteins = set(dedup_row["protein_in_syst"].split(","))
                if proteins & existing_proteins:
                    # subtypes match?
                    if dedup_row["subtype"] == row["subtype"]:
                        dedup_row["protein_in_syst"] = ",".join(sorted(existing_proteins | proteins))
                        dedup_row["tool_final"] = get_tool_final(
                            dedup_row["subtype"],
                            set([dedup_row["tool_final"], row["source"]]),
                            dedup_row["tool_final"]
                        )
                        overlap_found = True
                        successfully_deduplicated.append({
                            "seqid": genome,
                            "original_sys_id": row["sys_id"],
                            "merged_into_sys_id": dedup_row["sys_id"],
                            "subtype": dedup_row["subtype"],
                            "protein_in_syst": dedup_row["protein_in_syst"],
                            "tool_final": dedup_row["tool_final"]
                        })
                        break
                    else:
                        # check subtype resolution rule
                        resolved_subtype = resolve_subtype(dedup_row["subtype"], row["subtype"])
                        if resolved_subtype:
                            dedup_row["protein_in_syst"] = ",".join(sorted(existing_proteins | proteins))
                            dedup_row["tool_final"] = get_tool_final(
                                resolved_subtype,
                                set([dedup_row["tool_final"], row["source"]]),
                                dedup_row["tool_final"]
                            )
                            dedup_row["subtype"] = resolved_subtype
                            overlap_found = True
                            successfully_deduplicated.append({
                                "seqid": genome,
                                "original_sys_id": row["sys_id"],
                                "merged_into_sys_id": dedup_row["sys_id"],
                                "subtype": resolved_subtype,
                                "protein_in_syst": dedup_row["protein_in_syst"],
                                "tool_final": dedup_row["tool_final"]
                            })
                            break

        if not overlap_found:
            deduplicated.append({
                "seqid": row["seqid"],
                "sys_id": row["sys_id"],
                "type": row["type"],
                "subtype": row["subtype"],
                "activity": row["activity"],
                "sys_beg": row["sys_beg"],
                "sys_end": row["sys_end"],
                "protein_in_syst": row["protein_in_syst"],
                "genes_count": row["genes_count"],
                "name_of_profiles_in_sys": row["name_of_profiles_in_sys"],
                "tool_original": row["source"],
                "tool_final": row["source"]
            })

        # log any non-deduplicated overlaps with no rule
        for non_dedup_row in deduplicated:
            if non_dedup_row["seqid"] == genome and non_dedup_row["subtype"] != row["subtype"]:
                existing_proteins = set(non_dedup_row["protein_in_syst"].split(","))
                if proteins & existing_proteins:
                    if not resolve_subtype(non_dedup_row["subtype"], row["subtype"]):
                        non_deduplicated.append({
                            "seqid": genome,
                            "protein_overlap": ",".join(sorted(proteins & existing_proteins)),
                            "subtype_1": row["subtype"],
                            "subtype_2": non_dedup_row["subtype"]
                        })

# convert to DataFrames
deduplicated_df = pd.DataFrame(deduplicated)
successfully_deduplicated_df = pd.DataFrame(successfully_deduplicated)
non_deduplicated_df = pd.DataFrame(non_deduplicated)

# write outputs
deduplicated_df.to_csv(deduplicated_out, sep="\t", index=False)
successfully_deduplicated_df.to_csv(successfully_deduplicated_out, sep="\t", index=False)
non_deduplicated_df.to_csv(non_deduplicated_out, sep="\t", index=False)

print("✅ Deduplication complete:")
print(f"- deduplicated systems saved to deduplicated_defense_systems.tsv")
print(f"- non-deduplicated overlaps saved to non_deduplicated_overlaps.tsv")
print(f"- deduplication events saved to successfully_deduplicated_systems.tsv")

