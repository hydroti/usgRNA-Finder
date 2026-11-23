import pandas as pd
import re
import argparse
from dataclasses import dataclass, field
from typing import List
from enum import Enum
from openpyxl import Workbook
from openpyxl.styles import PatternFill, Font, Alignment
from openpyxl.utils import get_column_letter
from Bio import SeqIO


class ScissorMode(Enum):
    FORWARD = "Forward"
    REVERSE = "Reverse Complement"

@dataclass
class GeneInfo:

    GeneName: str
    Sequences: str
    originSequences: str

@dataclass
class GeneScissor:

    mode: ScissorMode
    pam_type: str
    geneName: str
    scissorString: str
    startIndex: int
    endIndex: int
    originStartIndex: int
    originEndIndex: int
    scissorScore: float = 0.0
    mmsIndexes: List[int] = field(default_factory=list)
    mms0Count: int = 0
    mms1Count: int = 0
    mms2Count: int = 0


PAM_PATTERNS = {
    "NGG": {
        "sense": lambda s: re.search(r"[ACGT]GG$", s),
        "antisense": lambda s: re.search(r"^CC[ACGT]", s)
    },
    "AGG": {
        "sense": lambda s: re.search(r"AGG$", s),
        "antisense": lambda s: re.search(r"^CCT", s)
    },
    "TGG": {
        "sense": lambda s: re.search(r"TGG$", s),
        "antisense": lambda s: re.search(r"^CCA", s)
    },
    "CGG": {
        "sense": lambda s: re.search(r"CGG$", s),
        "antisense": lambda s: re.search(r"^CCG", s)
    },
    "GGG": {
        "sense": lambda s: re.search(r"GGG$", s),
        "antisense": lambda s: re.search(r"^CCC", s)
    },
    "NG": {
        "sense": lambda s: re.search(r"[ACGT]G[ACGT]$", s),
        "antisense": lambda s: re.search(r"^[ACGT]C[ACGT]", s)
    },
    "AGN": {
        "sense": lambda s: re.search(r"AG[ATC]$", s),
        "antisense": lambda s: re.search(r"^[TAG]CT", s)
    },
    "TGN": {
        "sense": lambda s: re.search(r"TG[ATC]$", s),
        "antisense": lambda s: re.search(r"^[TAG]CA", s)
    },
    "GGN": {
        "sense": lambda s: re.search(r"GG[ATC]$", s),
        "antisense": lambda s: re.search(r"^[TAG]CC", s)
    },
    "CGN": {
        "sense": lambda s: re.search(r"CG[ATC]$", s),
        "antisense": lambda s: re.search(r"^[TAG]CG", s)
    },
    "AGA": {
        "sense": lambda s: re.search(r"AGA$", s),
        "antisense": lambda s: re.search(r"^TCT", s)
    },
    "AGT": {
        "sense": lambda s: re.search(r"AGT$", s),
        "antisense": lambda s: re.search(r"^ACT", s)
    },
    "AGC": {
        "sense": lambda s: re.search(r"AGC$", s),
        "antisense": lambda s: re.search(r"^GCT", s)
    },
    "TGA": {
        "sense": lambda s: re.search(r"TGA$", s),
        "antisense": lambda s: re.search(r"^TCA", s)
    },
    "TGT": {
        "sense": lambda s: re.search(r"TGT$", s),
        "antisense": lambda s: re.search(r"^ACA", s)
    },
    "TGC": {
        "sense": lambda s: re.search(r"TGC$", s),
        "antisense": lambda s: re.search(r"^GCA", s)
    },
    "GGA": {
        "sense": lambda s: re.search(r"GGA$", s),
        "antisense": lambda s: re.search(r"^TCC", s)
    },
    "GGT": {
        "sense": lambda s: re.search(r"GGT$", s),
        "antisense": lambda s: re.search(r"^ACC", s)
    },
    "GGC": {
        "sense": lambda s: re.search(r"GGC$", s),
        "antisense": lambda s: re.search(r"^GCC", s)
    },
    "CGA": {
        "sense": lambda s: re.search(r"CGA$", s),
        "antisense": lambda s: re.search(r"^TCG", s)
    },
    "CGT": {
        "sense": lambda s: re.search(r"CGT$", s),
        "antisense": lambda s: re.search(r"^ACG", s)
    },
    "CGC": {
        "sense": lambda s: re.search(r"CGC$", s),
        "antisense": lambda s: re.search(r"^GCG", s)
    }
}


CFD_MISMATCH_WEIGHTS = {

    ('A', 'A', 1): 1, ('A', 'A', 2): 0.727272727273, ('A', 'A', 3): 0.705882352941, ('A', 'A', 4): 0.636363636364,
    ('A', 'A', 5): 0.363636363636, ('A', 'A', 6): 0.714285714286, ('A', 'A', 7): 0.4375, ('A', 'A', 8): 0.428571428571,
    ('A', 'A', 9): 0.6, ('A', 'A', 10): 0.882352941176, ('A', 'A', 11): 0.307692307692, ('A', 'A', 12): 0.333333333333,
    ('A', 'A', 13): 0.3, ('A', 'A', 14): 0.533333333333, ('A', 'A', 15): 0.2, ('A', 'A', 16): 0,
    ('A', 'A', 17): 0.133333333333, ('A', 'A', 18): 0.5, ('A', 'A', 19): 0.538461538462, ('A', 'A', 20): 0.6,

    ('A', 'C', 1): 1, ('A', 'C', 2): 0.8, ('A', 'C', 3): 0.611111111111, ('A', 'C', 4): 0.625,
    ('A', 'C', 5): 0.72, ('A', 'C', 6): 0.714285714286, ('A', 'C', 7): 0.705882352941, ('A', 'C', 8): 0.733333333333,
    ('A', 'C', 9): 0.666666666667, ('A', 'C', 10): 0.555555555556, ('A', 'C', 11): 0.65, ('A', 'C', 12): 0.722222222222,
    ('A', 'C', 13): 0.652173913043, ('A', 'C', 14): 0.466666666667, ('A', 'C', 15): 0.65, ('A', 'C', 16): 0.192307692308,
    ('A', 'C', 17): 0.176470588235, ('A', 'C', 18): 0.4, ('A', 'C', 19): 0.375, ('A', 'C', 20): 0.764705882353,

    ('A', 'G', 1): 0.857142857143, ('A', 'G', 2): 0.785714285714, ('A', 'G', 3): 0.428571428571, ('A', 'G', 4): 0.352941176471,
    ('A', 'G', 5): 0.5, ('A', 'G', 6): 0.454545454545, ('A', 'G', 7): 0.4375, ('A', 'G', 8): 0.428571428571,
    ('A', 'G', 9): 0.571428571429, ('A', 'G', 10): 0.333333333333, ('A', 'G', 11): 0.4, ('A', 'G', 12): 0.263157894737,
    ('A', 'G', 13): 0.210526315789, ('A', 'G', 14): 0.214285714286, ('A', 'G', 15): 0.272727272727, ('A', 'G', 16): 0,
    ('A', 'G', 17): 0.176470588235, ('A', 'G', 18): 0.190476190476, ('A', 'G', 19): 0.206896551724, ('A', 'G', 20): 0.227272727273,

    ('C', 'A', 1): 1, ('C', 'A', 2): 0.909090909091, ('C', 'A', 3): 0.6875, ('C', 'A', 4): 0.8,
    ('C', 'A', 5): 0.636363636364, ('C', 'A', 6): 0.928571428571, ('C', 'A', 7): 0.8125, ('C', 'A', 8): 0.875,
    ('C', 'A', 9): 0.875, ('C', 'A', 10): 0.941176470588, ('C', 'A', 11): 0.307692307692, ('C', 'A', 12): 0.538461538462,
    ('C', 'A', 13): 0.7, ('C', 'A', 14): 0.733333333333, ('C', 'A', 15): 0.0666666666667, ('C', 'A', 16): 0.307692307692,
    ('C', 'A', 17): 0.466666666667, ('C', 'A', 18): 0.642857142857, ('C', 'A', 19): 0.461538461538, ('C', 'A', 20): 0.3,

    ('C', 'C', 1): 0.913043478261, ('C', 'C', 2): 0.695652173913, ('C', 'C', 3): 0.5, ('C', 'C', 4): 0.5,
    ('C', 'C', 5): 0.6, ('C', 'C', 6): 0.5, ('C', 'C', 7): 0.470588235294, ('C', 'C', 8): 0.642857142857,
    ('C', 'C', 9): 0.619047619048, ('C', 'C', 10): 0.388888888889, ('C', 'C', 11): 0.25, ('C', 'C', 12): 0.444444444444,
    ('C', 'C', 13): 0.136363636364, ('C', 'C', 14): 0, ('C', 'C', 15): 0.05, ('C', 'C', 16): 0.153846153846,
    ('C', 'C', 17): 0.0588235294118, ('C', 'C', 18): 0.133333333333, ('C', 'C', 19): 0.125, ('C', 'C', 20): 0.0588235294118,

    ('C', 'T', 1): 1, ('C', 'T', 2): 0.727272727273, ('C', 'T', 3): 0.866666666667, ('C', 'T', 4): 0.842105263158,
    ('C', 'T', 5): 0.571428571429, ('C', 'T', 6): 0.928571428571, ('C', 'T', 7): 0.75, ('C', 'T', 8): 0.65,
    ('C', 'T', 9): 0.857142857143, ('C', 'T', 10): 0.866666666667, ('C', 'T', 11): 0.75, ('C', 'T', 12): 0.714285714286,
    ('C', 'T', 13): 0.384615384615, ('C', 'T', 14): 0.35, ('C', 'T', 15): 0.222222222222, ('C', 'T', 16): 1,
    ('C', 'T', 17): 0.466666666667, ('C', 'T', 18): 0.538461538462, ('C', 'T', 19): 0.428571428571, ('C', 'T', 20): 0.5,

    ('G', 'A', 1): 1, ('G', 'A', 2): 0.636363636364, ('G', 'A', 3): 0.5, ('G', 'A', 4): 0.363636363636,
    ('G', 'A', 5): 0.3, ('G', 'A', 6): 0.666666666667, ('G', 'A', 7): 0.571428571429, ('G', 'A', 8): 0.625,
    ('G', 'A', 9): 0.533333333333, ('G', 'A', 10): 0.8125, ('G', 'A', 11): 0.384615384615, ('G', 'A', 12): 0.384615384615,
    ('G', 'A', 13): 0.3, ('G', 'A', 14): 0.266666666667, ('G', 'A', 15): 0.142857142857, ('G', 'A', 16): 0,
    ('G', 'A', 17): 0.25, ('G', 'A', 18): 0.666666666667, ('G', 'A', 19): 0.666666666667, ('G', 'A', 20): 0.7,

    ('G', 'G', 1): 0.714285714286, ('G', 'G', 2): 0.692307692308, ('G', 'G', 3): 0.384615384615, ('G', 'G', 4): 0.529411764706,
    ('G', 'G', 5): 0.785714285714, ('G', 'G', 6): 0.681818181818, ('G', 'G', 7): 0.6875, ('G', 'G', 8): 0.615384615385,
    ('G', 'G', 9): 0.538461538462, ('G', 'G', 10): 0.4, ('G', 'G', 11): 0.428571428571, ('G', 'G', 12): 0.529411764706,
    ('G', 'G', 13): 0.421052631579, ('G', 'G', 14): 0.428571428571, ('G', 'G', 15): 0.272727272727, ('G', 'G', 16): 0,
    ('G', 'G', 17): 0.235294117647, ('G', 'G', 18): 0.47619047619, ('G', 'G', 19): 0.448275862069, ('G', 'G', 20): 0.428571428571,

    ('G', 'T', 1): 0.9, ('G', 'T', 2): 0.846153846154, ('G', 'T', 3): 0.75, ('G', 'T', 4): 0.9,
    ('G', 'T', 5): 0.866666666667, ('G', 'T', 6): 1, ('G', 'T', 7): 1, ('G', 'T', 8): 1,
    ('G', 'T', 9): 0.642857142857, ('G', 'T', 10): 0.933333333333, ('G', 'T', 11): 1, ('G', 'T', 12): 0.933333333333,
    ('G', 'T', 13): 0.923076923077, ('G', 'T', 14): 0.75, ('G', 'T', 15): 0.941176470588, ('G', 'T', 16): 1,
    ('G', 'T', 17): 0.933333333333, ('G', 'T', 18): 0.692307692308, ('G', 'T', 19): 0.714285714286, ('G', 'T', 20): 0.9375,

    ('U', 'C', 1): 0.95652173913, ('U', 'C', 2): 0.84, ('U', 'C', 3): 0.5, ('U', 'C', 4): 0.625,
    ('U', 'C', 5): 0.64, ('U', 'C', 6): 0.571428571429, ('U', 'C', 7): 0.588235294118, ('U', 'C', 8): 0.733333333333,
    ('U', 'C', 9): 0.619047619048, ('U', 'C', 10): 0.5, ('U', 'C', 11): 0.4, ('U', 'C', 12): 0.5,
    ('U', 'C', 13): 0.260869565217, ('U', 'C', 14): 0, ('U', 'C', 15): 0.05, ('U', 'C', 16): 0.346153846154,
    ('U', 'C', 17): 0.117647058824, ('U', 'C', 18): 0.333333333333, ('U', 'C', 19): 0.25, ('U', 'C', 20): 0.176470588235,

    ('U', 'G', 1): 0.857142857143, ('U', 'G', 2): 0.857142857143, ('U', 'G', 3): 0.428571428571, ('U', 'G', 4): 0.647058823529,
    ('U', 'G', 5): 1, ('U', 'G', 6): 0.909090909091, ('U', 'G', 7): 0.6875, ('U', 'G', 8): 1,
    ('U', 'G', 9): 0.923076923077, ('U', 'G', 10): 0.533333333333, ('U', 'G', 11): 0.666666666667, ('U', 'G', 12): 0.947368421053,
    ('U', 'G', 13): 0.789473684211, ('U', 'G', 14): 0.285714285714, ('U', 'G', 15): 0.272727272727, ('U', 'G', 16): 0.666666666667,
    ('U', 'G', 17): 0.705882352941, ('U', 'G', 18): 0.428571428571, ('U', 'G', 19): 0.275862068966, ('U', 'G', 20): 0.0909090909091,

    ('U', 'T', 1): 1, ('U', 'T', 2): 0.846153846154, ('U', 'T', 3): 0.714285714286, ('U', 'T', 4): 0.47619047619,
    ('U', 'T', 5): 0.5, ('U', 'T', 6): 0.866666666667, ('U', 'T', 7): 0.875, ('U', 'T', 8): 0.8,
    ('U', 'T', 9): 0.928571428571, ('U', 'T', 10): 0.857142857143, ('U', 'T', 11): 0.75, ('U', 'T', 12): 0.8,
    ('U', 'T', 13): 0.692307692308, ('U', 'T', 14): 0.619047619048, ('U', 'T', 15): 0.578947368421, ('U', 'T', 16): 0.909090909091,
    ('U', 'T', 17): 0.533333333333, ('U', 'T', 18): 0.666666666667, ('U', 'T', 19): 0.285714285714, ('U', 'T', 20): 0.5625,
}


def read_excel_file(file_path: str) -> List[GeneInfo]:


    try:
        df = pd.read_excel(file_path, header=None, skiprows=1)
        gene_info_list = []
        for _, row in df.iterrows():
            if len(row) >= 2:
                sequences = str(row[1])
                origin_sequences = sequences.replace('-', '')
                gene_info = GeneInfo(
                    GeneName=str(row[0]),
                    Sequences=sequences,
                    originSequences=origin_sequences
                )
                gene_info_list.append(gene_info)
        return gene_info_list
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File '{file_path}' not found. Please ensure it is in the same directory as the script.")
    except Exception as e:
        raise Exception(f"Error occurred while reading the Excel file: {str(e)}")

def read_fasta_file(file_path):


    try:
        gene_info_list = []
        for record in SeqIO.parse(file_path, "fasta"):
            seq = str(record.seq).upper()
            origin_seq = seq.replace("-", "")
            gene_info_list.append(GeneInfo(
                GeneName=record.id,
                Sequences=seq,
                originSequences=origin_seq
            ))
        return gene_info_list
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File '{file_path}' not found. Please ensure it is in the same directory as the script.")
    except Exception as e:
        raise Exception(f"Error occurred while reading the FASTA file: {str(e)}")

def reverse_complement(dna_sequence: str) -> str:


    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-'}
    rc_sequence = "".join([complement_dict.get(base.upper(), base) for base in reversed(dna_sequence)])
    return rc_sequence

def contains_gap_or_ambiguous_base(sequence: str) -> bool:

    return bool(re.search(r"[^ACGT]", sequence.upper()))


def calc_cfd_mismatch_quality(guide: str, target: str, is_forward: bool) -> float:


    guide = guide.upper()
    target = target.upper()
    if len(guide) != 23 or len(target) != 23:
        raise ValueError("guide and target must both be 23 nt")
    if contains_gap_or_ambiguous_base(guide) or contains_gap_or_ambiguous_base(target):
        raise ValueError("guide and target must contain only A/C/G/T")

    if not is_forward:
        guide = reverse_complement(guide)
        target = reverse_complement(target)

    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    quality = 1.0
    for i in range(20):
        g = guide[i]
        t = target[i]
        if g != t:
            rna = 'U' if g == 'T' else g
            dna = complement_dict[t]
            weight = CFD_MISMATCH_WEIGHTS[(rna, dna, i + 1)]
            quality *= weight
    return quality


def find_scissors(target_gene: GeneInfo, pam_type: str) -> List[GeneScissor]:


    sequence = target_gene.Sequences.upper()
    scissors = []
    window_size = 23

    if pam_type not in PAM_PATTERNS:
        raise ValueError(f"Error: Unsupported PAM type '{pam_type}'. Supported types are: {', '.join(PAM_PATTERNS.keys())}")

    pam_patterns = PAM_PATTERNS[pam_type]

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        if contains_gap_or_ambiguous_base(window):
            continue


        if pam_patterns["sense"](window):
            origin_start = len(sequence[:i].replace('-', ''))
            origin_end = origin_start + len(window.replace('-', ''))
            scissors.append(GeneScissor(
                mode=ScissorMode.FORWARD,
                pam_type=pam_type,
                geneName=target_gene.GeneName,
                scissorString=window,
                startIndex=i,
                endIndex=i + window_size,
                originStartIndex=origin_start,
                originEndIndex=origin_end
            ))


        if pam_patterns["antisense"](window):
            origin_start = len(sequence[:i].replace('-', ''))
            origin_end = origin_start + len(window.replace('-', ''))
            scissors.append(GeneScissor(
                mode=ScissorMode.REVERSE,
                pam_type=pam_type,
                geneName=target_gene.GeneName,
                scissorString=window,
                startIndex=i,
                endIndex=i + window_size,
                originStartIndex=origin_start,
                originEndIndex=origin_end
            ))

    return scissors

def check_scissors_in_other_sequences(scissors: List[GeneScissor], other_genes: List[GeneInfo]) -> List[GeneScissor]:


    valid_scissors = []
    total_genes = len(other_genes)

    for scissor in scissors:
        scissor.mms0Count = 0
        scissor.mms1Count = 0
        scissor.mms2Count = 0
        scissor.mmsIndexes = []
        is_valid = True
        qualities = []

        pam_patterns = PAM_PATTERNS[scissor.pam_type]

        for gene in other_genes:

            if scissor.endIndex > len(gene.Sequences):
                is_valid = False
                break

            substring = gene.Sequences[scissor.startIndex:scissor.endIndex].upper()


            if len(substring) != 23 or contains_gap_or_ambiguous_base(substring):
                is_valid = False
                break


            if scissor.mode == ScissorMode.FORWARD:
                if not pam_patterns["sense"](substring):
                    is_valid = False
                    break
            else:
                if not pam_patterns["antisense"](substring):
                    is_valid = False
                    break


            mismatches = []
            is_forward = scissor.mode == ScissorMode.FORWARD


            for i in range(len(scissor.scissorString)):
                is_pam_position = (is_forward and i >= 20) or (not is_forward and i < 3)

                if not is_pam_position:
                    if scissor.scissorString[i].upper() != substring[i].upper():
                        mismatches.append(i + 1)


            if len(mismatches) > 2:
                is_valid = False
                break

            quality = calc_cfd_mismatch_quality(scissor.scissorString, substring, is_forward)
            qualities.append(quality)

            if len(mismatches) == 0:
                scissor.mms0Count += 1
            elif len(mismatches) == 1:
                scissor.mms1Count += 1
            elif len(mismatches) == 2:
                scissor.mms2Count += 1

            if mismatches:
                scissor.mmsIndexes.extend(mismatches)

        if not is_valid:
            continue

        scissor.scissorScore = round(sum(qualities) / total_genes * 100, 2)
        valid_scissors.append(scissor)

    return valid_scissors

def save_to_excel(scissors: List[GeneScissor], gene_info_list: List[GeneInfo], output_file: str):


    try:
        wb = Workbook()
        ws1 = wb.active
        ws1.title = "GeneScissor"
        ws2 = wb.create_sheet("GeneMatch")

        yellow_fill = PatternFill(start_color='FFFF00', end_color='FFFF00', fill_type='solid')

        headers1 = ["Rank", "PAM Type", "Position", "Strand", "Sequences(5'-3')", "Scores", "MM0", "MM1", "MM2"]
        for col, header in enumerate(headers1, 1):
            cell = ws1.cell(row=1, column=col, value=header)
            cell.font = Font(bold=True)
            cell.alignment = Alignment(horizontal='center')

        headers2 = ["Rank", "PAM Type", "Gene name", "Sequences(5'-3')", "MM Count"]
        for col, header in enumerate(headers2, 1):
            cell = ws2.cell(row=1, column=col, value=header)
            cell.font = Font(bold=True)
            cell.alignment = Alignment(horizontal='center')

        scissors.sort(key=lambda x: x.scissorScore, reverse=True)

        current_row = 2
        for rank, scissor in enumerate(scissors, 1):
            scissor_string = scissor.scissorString.upper()
            row_data = [
                rank,
                scissor.pam_type,
                f"{scissor.originStartIndex + 1}-{scissor.originEndIndex + 1}",
                '+' if scissor.mode == ScissorMode.FORWARD else '-',
                scissor_string,
                scissor.scissorScore,
                scissor.mms0Count,
                scissor.mms1Count,
                scissor.mms2Count
            ]

            for col, value in enumerate(row_data, 1):
                cell = ws1.cell(row=current_row, column=col, value=value)
                cell.fill = yellow_fill
            current_row += 1

        current_row = 2
        for rank, scissor in enumerate(scissors, 1):
            row_data = [rank, scissor.pam_type, scissor.geneName, scissor.scissorString.upper(), ""]
            for col, value in enumerate(row_data, 1):
                cell = ws2.cell(row=current_row, column=col, value=value)
                cell.fill = yellow_fill
            current_row += 1

            for gene in gene_info_list:
                if scissor.endIndex > len(gene.Sequences):
                    continue

                substring = gene.Sequences[scissor.startIndex:scissor.endIndex]

                if len(substring) != 23:
                    continue

                gene_string_list = list(substring)
                mismatch_count = 0

                is_forward = scissor.mode == ScissorMode.FORWARD


                if is_forward:
                    protospacer_substring = substring[:-3]
                else:
                    protospacer_substring = substring[3:]


                if '-' in protospacer_substring:
                    continue


                for i in range(len(gene_string_list)):
                    is_pam_position = (is_forward and i >= 20) or \
                                       (not is_forward and i < 3)

                    if not is_pam_position:
                        if gene_string_list[i].upper() != scissor.scissorString[i].upper():
                            gene_string_list[i] = gene_string_list[i].lower()
                            mismatch_count += 1
                        else:
                            gene_string_list[i] = gene_string_list[i].upper()
                    else:
                        gene_string_list[i] = gene_string_list[i].upper()

                gene_string = ''.join(gene_string_list)
                ws2.cell(row=current_row, column=1, value="")
                ws2.cell(row=current_row, column=2, value="")
                ws2.cell(row=current_row, column=3, value=gene.GeneName)
                ws2.cell(row=current_row, column=4, value=gene_string)
                ws2.cell(row=current_row, column=5, value=mismatch_count)
                current_row += 1

            current_row += 1

        for col in range(1, len(headers1) + 1):
            ws1.column_dimensions[get_column_letter(col)].width = 20
        for col in range(1, len(headers2) + 1):
            ws2.column_dimensions[get_column_letter(col)].width = 20

        wb.save(output_file)
    except Exception as e:
        raise Exception(f"Error occurred while saving the Excel file: {str(e)}")

def main():


    supported_pams = list(PAM_PATTERNS.keys())
    parser = argparse.ArgumentParser(description="usgRNA-Finder: A gene sequence analysis tool")
    parser.add_argument("--pam", type=str, default="NGG",
                        help=f"Specify the PAM type to search for. Supported types are: {', '.join(supported_pams)}. Default is NGG.")
    parser.add_argument('--input', type=str, default='input.xlsx', help='Input Excel or FASTA file')
    parser.add_argument('--output', type=str, default='output.xlsx', help='Output Excel file')
    args = parser.parse_args()

    try:
        if args.input.lower().endswith('.xlsx'):
            gene_info_list = read_excel_file(args.input)
        elif args.input.lower().endswith(('.fa', '.fasta')):
            gene_info_list = read_fasta_file(args.input)
        else:
            print("Error: Unsupported input file format. Please use a .xlsx or .fasta/.fa file.")
            return

        if not gene_info_list:
            print("Error: No gene data found.")
            return


        target_gene = gene_info_list[-1]
        print(f"\nFinding PAM type '{args.pam.upper()}' scissor sites from gene '{target_gene.GeneName}'...")
        scissors_all = find_scissors(target_gene, args.pam.upper())

        scissors_valid = check_scissors_in_other_sequences(scissors_all, gene_info_list)
        scissors_valid.sort(key=lambda x: x.scissorScore, reverse=True)

        print("\nValid scissor sites:")
        for scissor in scissors_valid:
            print(f"PAM Type: {scissor.pam_type}, Mode: {scissor.mode.value}")
            print(f"Gene: {scissor.geneName}")
            print(f"Sequence: {scissor.scissorString}")
            print(f"Position: {scissor.startIndex + 1}-{scissor.endIndex + 1}")
            print(f"Original Position: {scissor.originStartIndex + 1}-{scissor.originEndIndex + 1}")
            print(f"Score: {scissor.scissorScore:.2f}")
            mms_indices_str = ', '.join([str(x) for x in sorted(set(scissor.mmsIndexes))])
            print(f"Mismatch Positions: {mms_indices_str if mms_indices_str else 'None'}")
            print(f"Count with 0 mismatches: {scissor.mms0Count}")
            print(f"Count with 1 mismatch: {scissor.mms1Count}")
            print(f"Count with 2 mismatches: {scissor.mms2Count}")
            print("---")

        print(f"\nTotal scissor sites found: {len(scissors_all)}.")
        print(f"Valid scissor sites: {len(scissors_valid)}.")

        save_to_excel(scissors_valid, gene_info_list, args.output)
        print(f"\nResults saved to {args.output}.")

    except FileNotFoundError as e:
        print(e)
    except ValueError as e:
        print(e)
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    main()
