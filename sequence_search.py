import pandas as pd
import re
from dataclasses import dataclass, field
from typing import List, Set
from enum import Enum
from openpyxl import Workbook
from openpyxl.styles import PatternFill, Font, Alignment
from openpyxl.utils import get_column_letter

class ScissorMode(Enum):
    GG = "GG"  # GG ending mode
    CC = "CC"  # CC starting mode

@dataclass
class GeneInfo:
    GeneName: str  # Gene name
    Sequences: str  # Gene sequence
    originSequences: str  # Original gene sequence (sequence with '-' removed)

@dataclass
class GeneScissor:
    mode: ScissorMode  # GG ending or CC starting
    geneName: str  # Gene name
    scissorString: str  # Matched string
    startIndex: int  # Start position
    endIndex: int  # End position
    originStartIndex: int  # Start position in original sequence
    originEndIndex: int  # End position in original sequence
    scissorScore: int = 0  # Gene scissor score
    mmsIndexes: List[int] = field(default_factory=list)  # Mismatch positions
    mms0Count: int = 0  # Number of genes with 0 mismatches
    mms1Count: int = 0  # Number of genes with 1 mismatch
    mms2Count: int = 0  # Number of genes with 2 mismatches

def read_excel_file(file_path):
    # Read Excel file, skip header row
    df = pd.read_excel(file_path, header=None, skiprows=1)
    gene_info_list = []
    
    for _, row in df.iterrows():
        if len(row) >= 2:  # Ensure at least two columns of data
            sequences = str(row[1])
            # Remove '-' characters to get original sequence
            origin_sequences = sequences.replace('-', '')
            gene_info = GeneInfo(
                GeneName=str(row[0]),
                Sequences=sequences,
                originSequences=origin_sequences
            )
            gene_info_list.append(gene_info)
    
    return gene_info_list

def find_scissors_in_first_sequence(first_gene: GeneInfo) -> List[GeneScissor]:
    sequence = first_gene.Sequences
    origin_sequence = first_gene.originSequences
    scissors = []
    window_size = 23  # 21 letters + 2 GG/CC
    
    # Use sliding window to find patterns
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        
        # Check GG pattern
        if window.endswith('GG') and len(window[:-2]) == 21 and window[:-2].isalpha():
            # Calculate position in original sequence
            origin_start = len(sequence[:i].replace('-', ''))
            origin_end = origin_start + len(window.replace('-', ''))
            scissors.append(GeneScissor(
                mode=ScissorMode.GG,
                geneName=first_gene.GeneName,
                scissorString=window,
                startIndex=i,
                endIndex=i + window_size,
                originStartIndex=origin_start,
                originEndIndex=origin_end
            ))
        
        # Check CC pattern
        if window.startswith('CC') and len(window[2:]) == 21 and window[2:].isalpha():
            # Calculate position in original sequence
            origin_start = len(sequence[:i].replace('-', ''))
            origin_end = origin_start + len(window.replace('-', ''))
            scissors.append(GeneScissor(
                mode=ScissorMode.CC,
                geneName=first_gene.GeneName,
                scissorString=window,
                startIndex=i,
                endIndex=i + window_size,
                originStartIndex=origin_start,
                originEndIndex=origin_end
            ))
    
    return scissors

def calculate_position_score(mismatch_index: int, total_length: int, is_gg_pattern: bool) -> float:
    """Calculate score for mismatch position, higher score for closer to GG/CC"""
    if is_gg_pattern:
        distance = abs(mismatch_index - (total_length - 2))  # Distance to GG
    else:
        distance = abs(mismatch_index - 1)  # Distance to CC
    
    # Normalize score based on distance, closer is higher
    max_distance = total_length - 3  # Maximum possible distance
    normalized_score = 1 - (distance / max_distance)
    return normalized_score * 20  # Position score is 20 points

def check_scissors_in_other_sequences(scissors: List[GeneScissor], other_genes: List[GeneInfo]) -> List[GeneScissor]:
    valid_scissors = []
    total_genes = len(other_genes)
    
    for scissor in scissors:
        # Reset counters
        scissor.mms0Count = 0
        scissor.mms1Count = 0
        scissor.mms2Count = 0
        scissor.mmsIndexes = []
        
        is_valid = True
        skip_scissor = False
        
        # Determine if GG or CC pattern
        is_gg_pattern = scissor.mode == ScissorMode.GG
        
        for gene in other_genes:
            # Get substring at same position
            if scissor.endIndex > len(gene.Sequences):
                is_valid = False
                break
                
            substring = gene.Sequences[scissor.startIndex:scissor.endIndex]
            
            # Check if pattern matches
            if is_gg_pattern:
                if not re.match(r'[A-Za-z]{21}GG$', substring):
                    is_valid = False
                    break
            else:
                if not re.match(r'^CC[A-Za-z]{21}', substring):
                    is_valid = False
                    break
            
            # Check if other 20 letters contain "-"
            if is_gg_pattern:
                if '-' in substring[1:-2]:
                    skip_scissor = True
                    break
            else:
                if '-' in substring[2:-1]:
                    skip_scissor = True
                    break
            
            # Compare strings, record mismatch positions
            mismatches = []
            for i in range(len(scissor.scissorString)):
                # Skip first letter near CC/GG
                if (is_gg_pattern and i == len(substring)-3) or (not is_gg_pattern and i == 2):
                    continue
                if scissor.scissorString[i] != substring[i]:
                    mismatches.append(i)
            
            # Skip if more than 2 mismatches on any gene
            if len(mismatches) > 2:
                skip_scissor = True
                break
            
            # Record mismatch counts
            if len(mismatches) == 0:
                scissor.mms0Count += 1
            elif len(mismatches) == 1:
                scissor.mms1Count += 1
            elif len(mismatches) == 2:
                scissor.mms2Count += 1
            
            # Record mismatch positions
            scissor.mmsIndexes.extend(mismatches)
        
        if skip_scissor or not is_valid:
            continue
            
        # Calculate score (keep 2 decimal places)
        scissor.scissorScore = round(
            scissor.mms0Count / total_genes * 100 +
            scissor.mms1Count / total_genes * 80 +
            scissor.mms2Count / total_genes * 60,
            2
        )
        
        valid_scissors.append(scissor)
    
    return valid_scissors

def save_to_excel(scissors: List[GeneScissor], gene_info_list: List[GeneInfo], output_file: str):
    wb = Workbook()
    
    # Create two sheets
    ws1 = wb.active
    ws1.title = "基因剪刀"
    ws2 = wb.create_sheet("基因匹配")
    
    # Set styles
    yellow_fill = PatternFill(start_color='FFFF00', end_color='FFFF00', fill_type='solid')
    
    # Sheet1 headers
    headers1 = [
        "Rank", "Position", "Strand", "Sequences(5'-3')", "Scores", "MM0", "MM1", "MM2"
    ]
    
    # Set sheet1 header row
    for col, header in enumerate(headers1, 1):
        cell = ws1.cell(row=1, column=col, value=header)
        cell.font = Font(bold=True)
        cell.alignment = Alignment(horizontal='center')
    
    # Sheet2 headers
    headers2 = [
        "Rank", "Gene name", "Sequences(5'-3')", "MM Count"
    ]
    
    # Set sheet2 header row
    for col, header in enumerate(headers2, 1):
        cell = ws2.cell(row=1, column=col, value=header)
        cell.font = Font(bold=True)
        cell.alignment = Alignment(horizontal='center')
    
    # Sort scissors by score in descending order
    scissors.sort(key=lambda x: x.scissorScore, reverse=True)
    
    # Process sheet1
    current_row = 2
    for rank, scissor in enumerate(scissors, 1):
        # Process scissor string, all letters in uppercase
        scissor_string = scissor.scissorString.upper()
        
        # Write scissor info
        row_data = [
            rank,  # Rank (start from 1)
            f"{scissor.originStartIndex + 1}-{scissor.originEndIndex + 1}",  # Position (start from 1)
            '+' if scissor.mode == ScissorMode.GG else '-',  # Strand (+ for GG, - for CC)
            scissor_string,  # Sequences(5'-3') in uppercase
            scissor.scissorScore,  # Scores
            scissor.mms0Count,  # MM0
            scissor.mms1Count,  # MM1
            scissor.mms2Count   # MM2
        ]
        
        # Set row background color
        for col, value in enumerate(row_data, 1):
            cell = ws1.cell(row=current_row, column=col, value=value)
            cell.fill = yellow_fill
        
        current_row += 1
    
    # Process sheet2
    current_row = 2
    for rank, scissor in enumerate(scissors, 1):
        # Process scissor string, all letters in uppercase
        scissor_string = scissor.scissorString.upper()
        
        # Write scissor info
        row_data = [
            rank,  # Rank (start from 1)
            scissor.geneName,  # Gene name
            scissor_string,  # Sequences(5'-3') in uppercase
            ""  # MM Count (empty for scissor row)
        ]
        
        # Set row background color
        for col, value in enumerate(row_data, 1):
            cell = ws2.cell(row=current_row, column=col, value=value)
            cell.fill = yellow_fill
        
        current_row += 1
        
        # Write each gene's fragment
        for gene in gene_info_list:
            if scissor.endIndex > len(gene.Sequences):
                continue
                
            substring = gene.Sequences[scissor.startIndex:scissor.endIndex]
            # Process gene fragment, compare letters directly
            gene_string = list(substring)
            mismatch_count = 0
            for i in range(len(gene_string)):
                # Skip the first letter near GG/CC
                if (scissor.mode == ScissorMode.GG and i == len(gene_string)-3) or \
                   (scissor.mode == ScissorMode.CC and i == 2):
                    continue
                if gene_string[i].upper() == scissor.scissorString[i].upper():
                    gene_string[i] = gene_string[i].upper()  # Matching letters in uppercase
                else:
                    gene_string[i] = gene_string[i].lower()  # Mismatched letters in lowercase
                    mismatch_count += 1
            gene_string = ''.join(gene_string)
            
            # Write gene name, fragment and mismatch count
            ws2.cell(row=current_row, column=1, value="")  # Empty rank for gene fragments
            ws2.cell(row=current_row, column=2, value=gene.GeneName)
            ws2.cell(row=current_row, column=3, value=gene_string)
            ws2.cell(row=current_row, column=4, value=mismatch_count)
            
            current_row += 1
        
        # Add empty row
        current_row += 1
    
    # Adjust column width
    for col in range(1, len(headers1) + 1):
        ws1.column_dimensions[get_column_letter(col)].width = 15
    for col in range(1, len(headers2) + 1):
        ws2.column_dimensions[get_column_letter(col)].width = 15
    
    wb.save(output_file)

def main():
    try:
        # Read Excel file
        gene_info_list = read_excel_file('input.xlsx')
        
        if not gene_info_list:
            print("Error: No gene data found")
            return
        
        # Find scissor sites from the last gene sequence
        last_gene = gene_info_list[-1]
        
        print(f"\nFinding scissor sites from the last gene {last_gene.GeneName}...")
        scissors_all = find_scissors_in_first_sequence(last_gene)
        
        # Check scissor sites in all gene sequences
        scissors_valid = check_scissors_in_other_sequences(scissors_all, gene_info_list)
        
        # Sort by score
        scissors_valid.sort(key=lambda x: x.scissorScore, reverse=True)
        
        print("\nValid scissor sites:")
        for scissor in scissors_valid:
            print(f"Mode: {'+' if scissor.mode == ScissorMode.GG else '-'}")  # Use + for GG, - for CC
            print(f"Gene: {scissor.geneName}")
            print(f"Sequence: {scissor.scissorString}")
            print(f"Position: {scissor.startIndex + 1}-{scissor.endIndex + 1}")  # Start counting from 1
            print(f"Original Position: {scissor.originStartIndex + 1}-{scissor.originEndIndex + 1}")  # Start counting from 1
            print(f"Score: {scissor.scissorScore:.2f}")
            print(f"Mismatch Positions: {[x + 1 for x in sorted(set(scissor.mmsIndexes))]}")  # Start counting from 1
            print(f"0 Mismatch Count: {scissor.mms0Count}")
            print(f"1 Mismatch Count: {scissor.mms1Count}")
            print(f"2 Mismatch Count: {scissor.mms2Count}")
            print("---")
        
        print(f"\nTotal scissor sites found: {len(scissors_all)}")
        print(f"Valid scissor sites: {len(scissors_valid)}")
        
        # Save results to Excel file
        save_to_excel(scissors_valid, gene_info_list, 'output.xlsx')
        print("\nResults saved to output.xlsx")
        
    except FileNotFoundError:
        print("Error: input.xlsx file not found")
    except Exception as e:
        print(f"Error occurred: {str(e)}")

if __name__ == "__main__":
    main() 