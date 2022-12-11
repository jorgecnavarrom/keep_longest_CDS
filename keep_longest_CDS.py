#!/usr/bin/env python3

"""
Parses GenBank files and, if overlapping CDS features are detected (e.g. due to
isoforms), only keeps the longest one
"""

from pathlib import Path
from Bio import SeqIO
import argparse

__author__ = "Jorge Navarro"
__version__ = 1.1


def arg_parser():
    parser = argparse.ArgumentParser(description="CDS collapser. Removes \
        overlapping CDSs (if found) in GenBank files")

    parser.add_argument("-i", "--inputfolders", nargs='+', type=Path, \
        help="Folder(s) to search (recursively) for .gb, .gbk and .gbff files.")
    parser.add_argument("-f", "--files", nargs='+', type=Path, \
        help="Input individual files (accepted: .gb .gbk, .gbff.")

    default_output = Path("./") / "output"
    parser.add_argument("-o", "--outputfolder", default=default_output, 
        help=f"Base folder where results will be put " \
        f"(default='{default_output}').", type=Path)

    return parser.parse_args()


def get_file_list(individual_files: list, inputfolders: list) -> list:
    """Creates a single list of GenBank files
    """

    gbk_suffixes = {'gbk', 'gb', 'gbff'}
    all_found_files = set()

    try:
        for gbk in individual_files:
            if gbk.is_file() and gbk.suffix in gbk_suffixes:
                all_found_files.add(gbk)
    except TypeError:
        pass

    for folder in inputfolders:
        if folder.is_dir():
            for suffix in gbk_suffixes:
                for gbk in folder.rglob(f"*.{suffix}"):
                    all_found_files.add(gbk)

    return list(all_found_files)


def find_overlaps(CDS_list: list) -> list:
    """Finds groups of overlapping CDS

    CDS_list (list). Items: (start, end, CDS object)
    """

    regions = sorted(CDS_list, key=lambda tup: tup[0])
    group = [regions[0]]

    curr_start, curr_end, cds = regions[0]
    for region in regions[1:]:
        start, end, cds = region

        if start <= curr_end:
            curr_end = max(end, curr_end)
            group.append(region)
        else:
            yield group
            group = []

    if group:
        yield group


def de_overlap(file_list: list) -> tuple((dict, set)):
    """Removes overlapping CDS features to leave only the longest
    """
    
    modified = set()
    records = {}

    for gbk in file_list:
        try:
            gbk_recs = list(SeqIO.parse(str(gbk), "genbank"))
        except ValueError as e:
            print(f"Error, not able to parse  {str(gbk)}\n({str(e)})")
            continue

        for rec in gbk_recs:
            CDS_list = []

            for feature in rec.features:
                if feature.type != "CDS": continue
                
                CDS = feature
                
                cds_start = max(0, int(CDS.location.start))
                cds_end = max(0, int(CDS.location.end))

                CDS_list.append((cds_start, cds_end, CDS))

            # find all groups of overlapping CDS features in current record
            for group in find_overlaps(CDS_list):
                if len(group) < 2: continue

                # find longest CDS
                max_len = 0
                index_max = -1
                for idx, (start, end, CDS) in enumerate(group):
                    if len(CDS.qualifiers['translation'][0]) > max_len:
                        max_len = len(CDS.qualifiers['translation'][0])
                        index_max = idx

                # remove smaller CDSs
                for  idx, (start, end, CDS) in enumerate(group):
                    mark = ""
                    if idx == index_max: mark = ' <---'
                    print(f"{rec.id}\t{start}-{end} "\
                        f"{CDS.qualifiers['protein_id'][0]}\t"\
                        f"{len(CDS.qualifiers['translation'][0])} {mark}")

                    if idx != index_max:
                        # del CDS
                        rec.features.remove(CDS)
                        modified.add(gbk.stem)
                print()
        
        records[gbk.stem] = (gbk_recs)

    return records, modified


def write_records(outputfolder: Path, 
        processed_records: dict, 
        modified_records: set) -> None:

    outputfolder.mkdir(parents=True, exist_ok=True)

    for gbk_name, gbk_records in processed_records.items():
        modified_mark = ""
        if gbk_name in modified_records:
            modified_mark = '_mod'
        with open(outputfolder / f"{gbk_name}{modified_mark}.gbk", 'w') as f:
            SeqIO.write(gbk_records, f, "genbank")


def main():
    args = arg_parser()

    # collect files
    file_list = get_file_list(args.files, args.inputfolders)

    # process files and 
    processed_records, modified_records = de_overlap(file_list)

    # save genbank records in output
    write_records(args.outputfolder, processed_records, modified_records)


if __name__ == "__main__":
    main()