#!/usr/bin/env python3
"""
Convert gene-set GMT files from gene symbols to Entrez IDs

This script:
1. Detects if GMT files use symbols (alphanumeric) vs IDs (numeric)
2. Backs up symbol-based files with .original suffix
3. Converts symbols to Entrez IDs using gene reference file
4. Replaces original files with converted versions
5. Generates detailed conversion report
"""

import argparse
import os
import re
from pathlib import Path
from typing import Dict, List, Set, Tuple

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Convert GMT gene symbols to Entrez IDs")
    parser.add_argument("--gmt-files", required=True, help="Comma-separated list of GMT files")
    parser.add_argument("--gene-ref", required=True, help="Gene reference file (NCBI format)")
    parser.add_argument("--report", required=True, help="Output report file")
    parser.add_argument("--done", required=True, help="Output done marker file")
    return parser.parse_args()

def load_gene_mapping(gene_ref_file: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Load gene symbol to Entrez ID mapping from reference file
    
    Returns:
        symbol_to_id: Dict mapping gene symbols to Entrez IDs
        id_to_symbol: Dict mapping Entrez IDs to gene symbols
    """
    symbol_to_id = {}
    id_to_symbol = {}
    
    print(f"Loading gene reference from: {gene_ref_file}")
    
    with open(gene_ref_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                gene_id = parts[0].strip()
                gene_symbol = parts[5].strip()
                
                if gene_id and gene_symbol:
                    symbol_to_id[gene_symbol.upper()] = gene_id
                    id_to_symbol[gene_id] = gene_symbol
    
    print(f"Loaded {len(symbol_to_id)} gene symbol mappings")
    return symbol_to_id, id_to_symbol

def is_numeric_gene_id(gene: str) -> bool:
    """Check if gene identifier is numeric (Entrez ID)"""
    return gene.strip().isdigit()

def detect_gmt_type(gmt_file: str, sample_size: int = 10) -> str:
    """
    Detect if GMT file uses symbols or IDs by sampling genes
    
    Returns:
        'symbols' if genes are alphanumeric
        'entrez' if genes are numeric
        'mixed' if unclear
    """
    numeric_count = 0
    symbol_count = 0
    genes_checked = 0
    
    with open(gmt_file, 'r') as f:
        for line_num, line in enumerate(f):
            if line_num >= sample_size:
                break
            
            parts = line.strip().split('\t')
            if len(parts) > 2:
                # Check genes (skip name and description)
                genes = parts[2:]
                for gene in genes[:5]:  # Check first 5 genes per line
                    if gene.strip():
                        genes_checked += 1
                        if is_numeric_gene_id(gene):
                            numeric_count += 1
                        else:
                            symbol_count += 1
    
    if genes_checked == 0:
        return 'unknown'
    
    # If >80% are numeric, assume Entrez IDs
    if numeric_count / genes_checked > 0.8:
        return 'entrez'
    # If >80% are symbols, assume gene symbols
    elif symbol_count / genes_checked > 0.8:
        return 'symbols'
    else:
        return 'mixed'

def convert_gmt_file(gmt_file: str, symbol_to_id: Dict[str, str], 
                     report_lines: List[str]) -> Dict[str, any]:
    """
    Convert GMT file from symbols to Entrez IDs
    
    Returns dict with conversion statistics
    """
    stats = {
        'file': os.path.basename(gmt_file),
        'total_genesets': 0,
        'total_genes': 0,
        'converted_genes': 0,
        'unmapped_genes': 0,
        'skipped_genesets': 0,
        'unmapped_symbols': set()
    }
    
    # Check if already has backup
    backup_file = f"{gmt_file}.original"
    if os.path.exists(backup_file):
        print(f"  Backup already exists, skipping: {gmt_file}")
        report_lines.append(f"\n{stats['file']}: SKIPPED (backup exists)\n")
        stats['skipped'] = True
        return stats
    
    # Detect file type
    file_type = detect_gmt_type(gmt_file)
    
    if file_type == 'entrez':
        print(f"  Already using Entrez IDs, skipping: {gmt_file}")
        report_lines.append(f"\n{stats['file']}: SKIPPED (already Entrez IDs)\n")
        stats['skipped'] = True
        return stats
    
    if file_type == 'unknown':
        print(f"  WARNING: Cannot determine file type: {gmt_file}")
        report_lines.append(f"\n{stats['file']}: SKIPPED (unknown format)\n")
        stats['skipped'] = True
        return stats
    
    print(f"  Converting: {gmt_file}")
    stats['skipped'] = False
    
    # Read and convert
    converted_lines = []
    
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                converted_lines.append(line)
                continue
            
            stats['total_genesets'] += 1
            geneset_name = parts[0]
            description = parts[1]
            genes = parts[2:]
            
            converted_genes = []
            for gene in genes:
                gene = gene.strip()
                if not gene:
                    continue
                
                stats['total_genes'] += 1
                
                # If already numeric, keep it
                if is_numeric_gene_id(gene):
                    converted_genes.append(gene)
                    stats['converted_genes'] += 1
                else:
                    # Try to convert symbol to ID
                    gene_upper = gene.upper()
                    if gene_upper in symbol_to_id:
                        entrez_id = symbol_to_id[gene_upper]
                        converted_genes.append(entrez_id)
                        stats['converted_genes'] += 1
                    else:
                        # Unmapped - skip with warning
                        stats['unmapped_genes'] += 1
                        stats['unmapped_symbols'].add(gene)
            
            # Only keep gene-sets with at least one gene
            if converted_genes:
                converted_line = f"{geneset_name}\t{description}\t" + "\t".join(converted_genes)
                converted_lines.append(converted_line + "\n")
            else:
                stats['skipped_genesets'] += 1
    
    # Create backup
    print(f"  Creating backup: {backup_file}")
    os.rename(gmt_file, backup_file)
    
    # Write converted file
    with open(gmt_file, 'w') as f:
        f.writelines(converted_lines)
    
    # Add to report
    report_lines.append(f"\n{'='*80}\n")
    report_lines.append(f"File: {stats['file']}\n")
    report_lines.append(f"{'='*80}\n")
    report_lines.append(f"Total gene-sets: {stats['total_genesets']}\n")
    report_lines.append(f"Gene-sets kept: {stats['total_genesets'] - stats['skipped_genesets']}\n")
    report_lines.append(f"Gene-sets skipped (no valid genes): {stats['skipped_genesets']}\n")
    report_lines.append(f"Total genes processed: {stats['total_genes']}\n")
    report_lines.append(f"Genes successfully converted: {stats['converted_genes']}\n")
    report_lines.append(f"Genes unmapped (removed): {stats['unmapped_genes']}\n")
    
    if stats['unmapped_symbols']:
        report_lines.append(f"\nUnmapped symbols (sample, max 50):\n")
        for i, symbol in enumerate(sorted(stats['unmapped_symbols'])[:50]):
            report_lines.append(f"  - {symbol}\n")
        if len(stats['unmapped_symbols']) > 50:
            report_lines.append(f"  ... and {len(stats['unmapped_symbols']) - 50} more\n")
    
    return stats

def main():
    """Main execution"""
    args = parse_args()
    
    # Parse GMT file list
    gmt_files = [f.strip() for f in args.gmt_files.split(',')]
    
    print("\n" + "="*80)
    print("Gene-Set GMT Format Conversion")
    print("="*80)
    print(f"Processing {len(gmt_files)} GMT files")
    
    # Load gene mapping
    symbol_to_id, id_to_symbol = load_gene_mapping(args.gene_ref)
    
    # Process each file
    report_lines = []
    report_lines.append("="*80 + "\n")
    report_lines.append("Gene-Set Format Conversion Report\n")
    report_lines.append("="*80 + "\n")
    report_lines.append(f"\nTotal GMT files to process: {len(gmt_files)}\n")
    report_lines.append(f"Gene reference: {args.gene_ref}\n")
    report_lines.append(f"Total gene mappings available: {len(symbol_to_id)}\n")
    
    all_stats = []
    files_converted = 0
    files_skipped = 0
    
    for gmt_file in gmt_files:
        if not os.path.exists(gmt_file):
            print(f"  WARNING: File not found: {gmt_file}")
            report_lines.append(f"\n{os.path.basename(gmt_file)}: FILE NOT FOUND\n")
            continue
        
        stats = convert_gmt_file(gmt_file, symbol_to_id, report_lines)
        all_stats.append(stats)
        
        if stats.get('skipped', False):
            files_skipped += 1
        else:
            files_converted += 1
    
    # Summary
    total_genesets = sum(s['total_genesets'] for s in all_stats if not s.get('skipped', False))
    total_genes = sum(s['total_genes'] for s in all_stats if not s.get('skipped', False))
    total_converted = sum(s['converted_genes'] for s in all_stats if not s.get('skipped', False))
    total_unmapped = sum(s['unmapped_genes'] for s in all_stats if not s.get('skipped', False))
    
    report_lines.append(f"\n{'='*80}\n")
    report_lines.append(f"SUMMARY\n")
    report_lines.append(f"{'='*80}\n")
    report_lines.append(f"Files processed: {len(gmt_files)}\n")
    report_lines.append(f"Files converted: {files_converted}\n")
    report_lines.append(f"Files skipped: {files_skipped}\n")
    report_lines.append(f"Total gene-sets: {total_genesets}\n")
    report_lines.append(f"Total genes processed: {total_genes}\n")
    report_lines.append(f"Genes converted to Entrez IDs: {total_converted}\n")
    report_lines.append(f"Genes unmapped (removed): {total_unmapped}\n")
    
    if total_genes > 0:
        success_rate = (total_converted / total_genes) * 100
        report_lines.append(f"Conversion success rate: {success_rate:.1f}%\n")
    
    # Write report
    os.makedirs(os.path.dirname(args.report), exist_ok=True)
    with open(args.report, 'w') as f:
        f.writelines(report_lines)
    
    print("\n" + "="*80)
    print("Conversion Summary:")
    print(f"  Files converted: {files_converted}")
    print(f"  Files skipped: {files_skipped}")
    print(f"  Total genes processed: {total_genes}")
    print(f"  Genes converted: {total_converted}")
    print(f"  Genes unmapped: {total_unmapped}")
    if total_genes > 0:
        print(f"  Success rate: {success_rate:.1f}%")
    print(f"\nDetailed report: {args.report}")
    print("="*80 + "\n")
    
    # Create done marker
    os.makedirs(os.path.dirname(args.done), exist_ok=True)
    with open(args.done, 'w') as f:
        f.write(f"Format genesets completed\n")
        f.write(f"Files converted: {files_converted}\n")
        f.write(f"Files skipped: {files_skipped}\n")
    
    print(f"Done marker created: {args.done}")

if __name__ == "__main__":
    main()
