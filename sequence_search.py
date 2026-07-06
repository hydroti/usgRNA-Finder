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

# Enum for PAM mode and strand direction
class ScissorMode(Enum):
    FORWARD = "Forward"
    REVERSE = "Reverse Complement"

@dataclass
class GeneInfo:
    """Dataclass to store gene name and sequence information"""
    GeneName: str
    Sequences: str
    originSequences: str

@dataclass
class GeneScissor:
    """Dataclass to store gene scissor site information"""
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

# Define all supported PAM types and their forward and reverse complement lambda functions
# Each function takes a 23-base sequence string as input.
# 'sense' function checks for the PAM pattern at the 3' end.
# 'antisense' function checks for the PAM pattern at the 5' end.
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
        "sense": lambda s: re.search(r"[ACGT]G[ATC]$", s),
        "antisense": lambda s: re.search(r"^[TAG]G[ACGT]", s)
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

def read_excel_file(file_path: str) -> List[GeneInfo]:
    """
    Reads an Excel file, skips the header row, and extracts gene names and sequences.
    """
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
    """
    Reads a FASTA file and extracts gene names and sequences.
    """
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
    """
    Calculates the reverse complement of a DNA sequence.
    """
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', '-': '-'}
    rc_sequence = "".join([complement_dict.get(base.upper(), base) for base in reversed(dna_sequence)])
    return rc_sequence

def find_scissors(target_gene: GeneInfo, pam_type: str) -> List[GeneScissor]:
    """
    Finds gene scissor sites in the target gene sequence based on the specified PAM type.
    """
    sequence = target_gene.Sequences.upper()
    scissors = []
    window_size = 23
    
    if pam_type not in PAM_PATTERNS:
        raise ValueError(f"Error: Unsupported PAM type '{pam_type}'. Supported types are: {', '.join(PAM_PATTERNS.keys())}")
    
    pam_patterns = PAM_PATTERNS[pam_type]

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        if '-' in window:
            continue
        
        # Search for forward strand (NGG-like)
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

        # Search for reverse complement strand (CCN-like) on the same input sequence
        if pam_patterns["antisense"](window):
            origin_start = len(sequence[:i].replace('-', ''))
            origin_end = origin_start + len(window.replace('-', ''))
            scissors.append(GeneScissor(
                mode=ScissorMode.REVERSE,
                pam_type=pam_type,
                geneName=target_gene.GeneName,
                scissorString=window,  # Store the raw sequence from the input strand
                startIndex=i,
                endIndex=i + window_size,
                originStartIndex=origin_start,
                originEndIndex=origin_end
            ))
            
    return scissors

def check_scissors_in_other_sequences(scissors: List[GeneScissor], other_genes: List[GeneInfo]) -> List[GeneScissor]:
    """
    Checks the validity of gene scissor sites in all other gene sequences and calculates scores.
    This version has been updated to include a strict PAM pattern and gap check in the protospacer.
    """ 
    valid_scissors = []
    total_genes = len(other_genes)

    for scissor in scissors:
        scissor.mms0Count = 0
        scissor.mms1Count = 0
        scissor.mms2Count = 0
        scissor.mmsIndexes = []
        is_valid = True

        pam_patterns = PAM_PATTERNS[scissor.pam_type]

        for gene in other_genes:
            # 1. Check if the gene sequence is long enough to contain the scissor site
            if scissor.endIndex > len(gene.Sequences):
                is_valid = False
                break
            
            substring = gene.Sequences[scissor.startIndex:scissor.endIndex].upper()
            
            # 2. Check if the substring has the correct length
            if len(substring) != 23:
                is_valid = False
                break
            
            # --- START: Modified Logic ---
            # 3. Re-validate the PAM pattern
            if scissor.mode == ScissorMode.FORWARD:
                if not pam_patterns["sense"](substring):
                    is_valid = False
                    break
            else: # Reverse Complement
                if not pam_patterns["antisense"](substring):
                    is_valid = False
                    break

            # 4. Check the protospacer region for gaps ('-')
            if scissor.mode == ScissorMode.FORWARD:
                protospacer_substring = substring[:-3]
            else:
                protospacer_substring = substring[3:]
            
            if '-' in protospacer_substring:
                is_valid = False
                break
            # --- END: Modified Logic ---
            
            mismatches = []
            is_forward = scissor.mode == ScissorMode.FORWARD
            
            # Compare sequences to count mismatches, skipping the PAM region
            for i in range(len(scissor.scissorString)):
                is_pam_position = (is_forward and i >= 20) or (not is_forward and i < 3)
                
                if not is_pam_position:
                    if scissor.scissorString[i].upper() != substring[i].upper():
                        mismatches.append(i + 1)
            
            # Skip if more than 2 mismatches on any gene
            if len(mismatches) > 2:
                is_valid = False
                break
            
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

        scissor.scissorScore = round(
            scissor.mms0Count / total_genes * 100 +
            scissor.mms1Count / total_genes * 80 +
            scissor.mms2Count / total_genes * 60,
            2
        )
        valid_scissors.append(scissor)
        
    return valid_scissors

def save_to_excel(scissors: List[GeneScissor], gene_info_list: List[GeneInfo], output_file: str):
    """
    Saves the analysis results to an Excel file with two sheets.
    """
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
                
                # Check protospacer for '-'
                if is_forward:
                    protospacer_substring = substring[:-3]
                else:
                    protospacer_substring = substring[3:]

                # If protospacer contains '-', skip to the next gene
                if '-' in protospacer_substring:
                    continue
                
                # Mark mismatches and count them
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
    """
    Main function to handle command-line arguments and run the analysis.
    """
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

        # Find scissor sites from the last gene sequence
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

