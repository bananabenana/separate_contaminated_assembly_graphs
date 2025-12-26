#!/usr/bin/env python3

import os
import glob
import re
import argparse
import subprocess
from collections import defaultdict, Counter, defaultdict
import networkx as nx

def parse_args():
    parser = argparse.ArgumentParser(description="Split assemblies into membership groups by taxon.")
    parser.add_argument("--gfa_directory", required=True, help="Directory containing GFA files")
    parser.add_argument("--fasta_dir", required=True, help="Directory containing FASTA files")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--taxa_level", required=True, choices=["kingdom","phylum","class","order","family","genus","species"])
    parser.add_argument("--threads", default="2", help="Threads for Kraken2")
    parser.add_argument("--kraken2_db", required=True, help="Kraken2 DB path")
    parser.add_argument("--taxonkit_db", required=True, help="Taxonkit DB path")
    return parser.parse_args()


def natural_key_sorting(text):
    # split into digit and non-digit chunks: 'contig10a' -> ['contig', 10, 'a']
    return [int(t) if t.isdigit() else t.lower()
            for t in re.split(r'(\d+)', text)]


def read_gfa_links(gfa_file):
    contigs = set()
    links = []
    with open(gfa_file) as f:
        for line in f:
            if line.startswith("S"):
                contigs.add(line.strip().split("\t")[1])
            elif line.startswith("L"):
                parts = line.strip().split("\t")
                links.append((parts[1], parts[3]))
    return contigs, links


def connected_components(contigs, links):
    G = nx.Graph()
    G.add_nodes_from(contigs)
    G.add_edges_from(links)
    return [set(c) for c in nx.connected_components(G)]


def append_basename_to_fasta_and_gfa(fasta_file, gfa_file, suffix, fasta_out, gfa_out):
    # Append to FASTA
    with open(fasta_file) as f, open(fasta_out, "w") as outf:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].strip().split(maxsplit=1)  # split at first whitespace
                contig_id = parts[0] + f"___{suffix}"
                rest = " " + parts[1] if len(parts) > 1 else ""
                outf.write(f">{contig_id}{rest}\n")
            else:
                outf.write(line)

    # Append to GFA
    with open(gfa_file) as f, open(gfa_out, "w") as outg:
        for line in f:
            if line.startswith("S"):
                parts = line.strip().split("\t")
                parts[1] = parts[1] + f"___{suffix}"
                outg.write("\t".join(parts) + "\n")
            elif line.startswith("L"):
                parts = line.strip().split("\t")
                parts[1] = parts[1] + f"___{suffix}"
                parts[3] = parts[3] + f"___{suffix}"
                outg.write("\t".join(parts) + "\n")
            else:
                outg.write(line)


def run_kraken2(fasta_files, kraken2_db, threads, output_file):
    cmd = ["kraken2", "--threads", str(threads), "--db", kraken2_db] + fasta_files
    with open(output_file, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)


def prep_taxonkit_input(kraken_output, taxonkit_input):
    with open(kraken_output) as inp, open(taxonkit_input, "w") as out:
        for line in inp:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            contig = cols[1]
            taxid = cols[2]
            out.write(f"{contig}\t{taxid}\n")


def run_taxonkit_lineage(input_file, taxonkit_db, output_file):
    cmd = ["taxonkit", "lineage", input_file, "-i", "2", "--delimiter", "	", "--data-dir", taxonkit_db]
    with open(output_file, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)


def read_taxonkit_lineage(lineage_file):
    # contig_lineage: contig_name -> dict of ranks
    contig_lineage = {}
    with open(lineage_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            taxid = cols[1]
            contig = None  # we need to recover this from the prep mapping
            # Then assign lineage
            lineage_dict = {
                "kingdom": cols[3] if len(cols)>3 else "",
                "phylum": cols[4] if len(cols)>4 else "",
                "class": cols[5] if len(cols)>5 else "",
                "order": cols[6] if len(cols)>6 else "",
                "family": cols[7] if len(cols)>7 else "",
                "genus": cols[9] if len(cols)>9 else "",
                "species": cols[11] if len(cols)>11 else "",
            }
            contig_lineage[taxid] = lineage_dict
    return contig_lineage


def merge_memberships_by_taxa(membership_dict, contig_lineage, taxa_level):
    """
    membership_dict: dict of membership_id -> set of contigs
    contig_lineage: dict contig -> lineage dict
    taxa_level: rank to use for grouping
    Returns: dict of new_group_id -> set of contigs
    """
    taxa_to_memberships = defaultdict(list)
    # assign each membership to a taxon at the selected level
    membership_taxon = {}
    for mem_id, contigs in membership_dict.items():
        taxa = [contig_lineage[c].get(taxa_level,"NA") for c in contigs if c in contig_lineage]
        if taxa:
            # most common taxon in membership
            taxon = Counter(taxa).most_common(1)[0][0]
        else:
            taxon = "NA"
        membership_taxon[mem_id] = taxon
        taxa_to_memberships[taxon].append(mem_id)

    # merge memberships with same taxon
    merged = {}
    new_id = 1
    for taxon, mem_ids in taxa_to_memberships.items():
        merged_contigs = set()
        for mem_id in mem_ids:
            merged_contigs.update(membership_dict[mem_id])
        merged[new_id] = merged_contigs
        new_id += 1
    return merged


def write_membership_outputs(base_name, components, contig_lineage,
                             appended_fasta_file, appended_gfa_file, output_dir):
    """
    Write FASTA/GFA/TSV per membership:
        FASTA: only sequences in the membership
        GFA:   only S lines for those contigs + L lines whose two endpoints are both in the membership
        TSV:   with contig_length column, sorted by contig
    """
    import os

    print(f"[INFO] Writing membership outputs for {base_name}")

    # Load fasta sequences and lengths
    fasta_seqs, seq_lengths = {}, {}
    with open(appended_fasta_file) as f:
        current, seq_lines = None, []
        for line in f:
            if line.startswith(">"):
                if current:
                    seq = "".join(seq_lines)
                    fasta_seqs[current] = seq
                    seq_lengths[current] = len(seq)
                current = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if current:
            seq = "".join(seq_lines)
            fasta_seqs[current] = seq
            seq_lengths[current] = len(seq)

    # Index GFA lines just once
    S_lines = {}           # contig -> S line string
    L_by_node = {}         # node -> set of link-line strings it participates in
    with open(appended_gfa_file) as f:
        for line in f:
            if line.startswith("S"):
                parts = line.rstrip("\n").split("\t")
                S_lines[parts[1]] = parts
            elif line.startswith("L"):
                p = line.rstrip("\n").split("\t")
                L_by_node.setdefault(p[1], set()).add(tuple(p))
                L_by_node.setdefault(p[3], set()).add(tuple(p))

    # Write per-membership outputs
    for i, contig_set in components.items():
        print(f"[INFO]   Membership {i}: {len(contig_set)} contigs")
        comp_fasta = os.path.join(output_dir, f"{base_name}_membership_{i}.fasta")
        comp_gfa   = os.path.join(output_dir, f"{base_name}_membership_{i}.gfa")
        comp_tsv   = os.path.join(output_dir, f"{base_name}_membership_{i}.tsv")

        rows = []

        # FASTA
        with open(comp_fasta, "w") as outf:
            for c in contig_set:
                if c in fasta_seqs:
                    clean = c.split("___")[0]
                    outf.write(f">{clean}\n{fasta_seqs[c]}\n")

        # GFA: only once per S or L line
        with open(comp_gfa, "w") as outg:
            # S lines
            for c in contig_set:
                if c in S_lines:
                    parts = S_lines[c].copy()
                    parts[1] = parts[1].split("___")[0]
                    outg.write("\t".join(parts) + "\n")
            # L lines where both ends are in contig_set
            seen_links = set()
            for c in contig_set:
                for p in L_by_node.get(c, []):
                    if p[1] in contig_set and p[3] in contig_set:
                        if p not in seen_links:
                            seen_links.add(p)
                            q = list(p)
                            q[1] = q[1].split("___")[0]
                            q[3] = q[3].split("___")[0]
                            outg.write("\t".join(q) + "\n")

        # TSV rows
        for c in contig_set:
            clean = c.split("___")[0]
            lin = contig_lineage.get(c, {})
            rows.append([
                base_name,
                clean,
                str(i),
                str(seq_lengths.get(c, "")),
                lin.get("kingdom",""), lin.get("phylum",""), lin.get("class",""),
                lin.get("order",""),   lin.get("family",""),
                lin.get("genus",""),   lin.get("species","")
            ])
        rows.sort(key=lambda r: natural_key_sorting(r[1]))
        with open(comp_tsv, "w") as tsv:
            tsv.write("\t".join([
                "input_filename","contig","membership","contig_length",
                "kingdom","phylum","class","order","family","genus","species"
            ]) + "\n")
            for r in rows:
                tsv.write("\t".join(r) + "\n")


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    gfa_files = glob.glob(os.path.join(args.gfa_directory, "*.gfa"))
    appended_fastas = []

    # Append basename to contigs for Kraken2
    temp_files = []
    for gfa_file in gfa_files:
        base = os.path.splitext(os.path.basename(gfa_file))[0]
        fasta_file = os.path.join(args.fasta_dir, f"{base}.fasta")
        if not os.path.exists(fasta_file):
            print(f"Skipping {base}, fasta not found")
            continue
        temp_fasta = os.path.join(args.output_dir, f"{base}_appended.fasta")
        temp_gfa = os.path.join(args.output_dir, f"{base}_appended.gfa")
        append_basename_to_fasta_and_gfa(fasta_file, gfa_file, base, temp_fasta, temp_gfa)
        appended_fastas.append(temp_fasta)
        temp_files.append((base, temp_fasta, temp_gfa, fasta_file))  # store info for later

    # Run Kraken2 once on all appended FASTAs
    print("[INFO] Running Kraken2 ...")
    kraken_output = os.path.join(args.output_dir, "all_kraken_output.tsv")
    run_kraken2(appended_fastas, args.kraken2_db, args.threads, kraken_output)

    # Prepare Taxonkit input (contig -> taxid)
    print("[INFO] Preparing Taxonkit input ...")
    taxonkit_input = os.path.join(args.output_dir, "taxonkit_prep.tsv")
    prep_taxonkit_input(kraken_output, taxonkit_input)

    # Run Taxonkit lineage
    print("[INFO] Running Taxonkit lineage command ...")
    lineage_file = os.path.join(args.output_dir, "taxonkit_lineage.tsv")
    run_taxonkit_lineage(taxonkit_input, args.taxonkit_db, lineage_file)

    # Parse contig -> taxid mapping
    print("[INFO] Parsing contig to taxid mapping ...")
    contig_taxid_map = {}
    with open(taxonkit_input) as f:
        for line in f:
            contig, taxid = line.strip().split("\t")
            contig_taxid_map[contig] = taxid

    # Parse Taxonkit output (taxid -> full lineage)
    print(f"[INFO] Parsing taxid to lineage mapping ...")
    taxid_lineage_map = {}
    with open(lineage_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            taxid = cols[1]
            lineage_dict = {
                "kingdom": cols[3] if len(cols) > 3 else "",
                "phylum": cols[4] if len(cols) > 4 else "",
                "class": cols[5] if len(cols) > 5 else "",
                "order": cols[6] if len(cols) > 6 else "",
                "family": cols[7] if len(cols) > 7 else "",
                "genus": cols[9] if len(cols) > 9 else "",
                "species": cols[11] if len(cols) > 11 else "",
            }
            taxid_lineage_map[taxid] = lineage_dict

    # Build contig -> full lineage mapping
    print(f"[INFO] Building contig to full lineage mapping ...")
    contig_lineage = {}
    for contig, taxid in contig_taxid_map.items():
        contig_lineage[contig] = taxid_lineage_map.get(taxid, {})

    
    # Merge memberships and write outputs
    print(f"[INFO] Merging contig membership groups ...")
    for base, temp_fasta, temp_gfa, orig_fasta in temp_files:
        contigs, links = read_gfa_links(temp_gfa)
        components = connected_components(contigs, links)
        membership_dict = {i + 1: comp for i, comp in enumerate(components)}
        merged_memberships = merge_memberships_by_taxa(
            membership_dict, contig_lineage, args.taxa_level
        )
        print(f"[INFO] {base}: {len(membership_dict)} initial memberships, "
              f"{len(merged_memberships)} after merging by {args.taxa_level}")

        write_membership_outputs(
            base,
            merged_memberships,
            contig_lineage,
            temp_fasta,
            temp_gfa,
            args.output_dir
        )

        # Remove temporary appended files
        try:
            os.remove(temp_fasta)
            os.remove(temp_gfa)
            print(f"[INFO] Cleaned temporary files for {base}")
        except OSError as e:
            print(f"[WARN] Could not remove temp files for {base}: {e}")

    print("Done!")

if __name__=="__main__":
    main()
