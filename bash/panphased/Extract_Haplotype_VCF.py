#!/usr/bin/env python3
"""
Extract a specific haplotype from a sample in a phased VCF file.

Usage:
    python extract_haplo.py SAMPLE HAPLOTYPE INPUT_VCF [OUTPUT_VCF]

Arguments:
    SAMPLE      - Sample ID to extract from
    HAPLOTYPE   - Haplotype number to extract (1-based)
    INPUT_VCF   - Input VCF file (can be gzipped)
    OUTPUT_VCF  - Output VCF file (optional, default: SAMPLE_HAPLOTYPE.vcf.gz)

Options:
    --keep-format  Keep all FORMAT fields (default: only keep GT field)
    --keep-info    Keep INFO field data (default: replace with ".")

Example:
    python extract_haplo.py H6 2 input.vcf.gz
    This will extract haplotype 2 from sample H6 and save it as H6_2.vcf.gz
"""

import sys
import gzip
import os
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Extract a specific haplotype from a sample in a phased VCF file.')
    parser.add_argument('sample', help='Sample ID to extract from')
    parser.add_argument('haplotype', type=int, help='Haplotype number to extract (1-based)')
    parser.add_argument('input_vcf', help='Input VCF file (can be gzipped)')
    parser.add_argument('output_vcf', nargs='?', default=None, 
                        help='Output VCF file (default: SAMPLE_HAPLOTYPE.vcf.gz)')
    parser.add_argument('--keep-format', action='store_true', 
                        help='Keep all FORMAT fields (default: only keep GT field)')
    parser.add_argument('--keep-info', action='store_true',
                        help='Keep INFO field data (default: replace with ".")')
    
    args = parser.parse_args()
    
    # Validate haplotype number is positive
    if args.haplotype < 1:
        parser.error("Haplotype number must be a positive integer (1-based)")
    
    # Set default output filename if not provided
    if args.output_vcf is None:
        base = f"{args.sample}_{args.haplotype}"
        args.output_vcf = f"{base}.vcf.gz"
        
    return args

def is_gzipped(filename):
    """Check if a file is gzipped based on its extension."""
    return filename.endswith('.gz')

def open_vcf(filename, mode='r'):
    """Open a VCF file, handling both gzipped and plain text files."""
    if is_gzipped(filename):
        if 'b' not in mode:
            mode += 'b'
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def process_vcf(input_vcf, output_vcf, sample_id, haplotype_idx, keep_format=False, keep_info=False):
    """
    Process the VCF file to extract a specific haplotype from a sample.
    
    Args:
        input_vcf: Path to input VCF file
        output_vcf: Path to output VCF file
        sample_id: ID of the sample to extract from
        haplotype_idx: 1-based index of the haplotype to extract
        keep_format: If True, keep all FORMAT fields; if False, only keep GT
        keep_info: If True, keep INFO field data; if False, replace with "."
    """
    # Convert to 0-based index for internal use
    haplotype_idx = haplotype_idx - 1
    
    # Determine if output should be gzipped
    output_is_gzipped = is_gzipped(output_vcf)
    
    # Open input file
    with open_vcf(input_vcf) as in_file:
        # Open output file with appropriate mode
        if output_is_gzipped:
            out_file = gzip.open(output_vcf, 'wb')
        else:
            out_file = open(output_vcf, 'w')
        
        try:
            sample_found = False
            sample_col_idx = -1
            
            # Process each line
            for line in in_file:
                if isinstance(line, bytes):
                    line = line.decode('utf-8')
                
                # Handle header lines
                if line.startswith('#'):
                    if line.startswith('#CHROM'):
                        # Find the sample column
                        header_fields = line.strip().split('\t')
                        if sample_id in header_fields:
                            sample_col_idx = header_fields.index(sample_id)
                            sample_found = True
                            
                            # Create new header with just the target sample (renamed)
                            new_header = header_fields[:9]  # Keep standard VCF columns
                            new_header.append(f"{sample_id}_{haplotype_idx+1}")  # Add renamed sample
                            new_header_line = '\t'.join(new_header) + '\n'
                            
                            # Write in the appropriate format
                            if output_is_gzipped:
                                out_file.write(new_header_line.encode('utf-8'))
                            else:
                                out_file.write(new_header_line)
                        else:
                            sys.stderr.write(f"Error: Sample '{sample_id}' not found in VCF header.\n")
                            sys.stderr.write(f"Available samples: {', '.join(header_fields[9:])}\n")
                            sys.exit(1)
                    elif line.startswith('##FORMAT=') and not keep_format:
                        # If we're not keeping all FORMAT fields, only keep the GT field
                        if line.startswith('##FORMAT=<ID=GT,'):
                            if output_is_gzipped:
                                out_file.write(line.encode('utf-8'))
                            else:
                                out_file.write(line)
                        # Skip other FORMAT fields
                    elif line.startswith('##INFO=') and not keep_info:
                        # Skip INFO headers if we're not keeping INFO
                        continue
                    else:
                        # Pass through other header lines unchanged
                        if output_is_gzipped:
                            if isinstance(line, str):
                                out_file.write(line.encode('utf-8'))
                            else:
                                out_file.write(line)
                        else:
                            if isinstance(line, bytes):
                                out_file.write(line.decode('utf-8'))
                            else:
                                out_file.write(line)
                else:
                    # Process variant lines
                    if not sample_found:
                        sys.stderr.write("Error: No sample column found in VCF header.\n")
                        sys.exit(1)
                    
                    fields = line.strip().split('\t')
                    
                    if sample_col_idx >= len(fields):
                        sys.stderr.write(f"Error: Invalid VCF format. Line has fewer columns than expected.\n")
                        continue
                    
                    # Extract genotype data for the sample
                    gt_data = fields[sample_col_idx]
                    gt_fields = gt_data.split(':')
                    format_fields = fields[8].split(':')
                    
                    if not gt_fields:
                        sys.stderr.write(f"Warning: Empty genotype data for sample at line: {line.strip()}\n")
                        continue
                    
                    # Extract the GT field (should be first)
                    gt = gt_fields[0]
                    
                    # For phased data, split by pipe
                    if '|' in gt:
                        alleles = gt.split('|')
                        
                        # Check if haplotype index is valid
                        if haplotype_idx >= len(alleles):
                            sys.stderr.write(f"Warning: Haplotype {haplotype_idx+1} does not exist for sample {sample_id} at position {fields[0]}:{fields[1]}. Ploidy is {len(alleles)}. Skipping.\n")
                            continue
                        
                        # Create new genotype with single haplotype (haploid format)
                        selected_allele = alleles[haplotype_idx]
                        new_gt = selected_allele  # Just use the selected allele as-is for haploid output
                    else:
                        # Handle unphased data (though this shouldn't happen in phased VCF)
                        sys.stderr.write(f"Warning: Unphased genotype found at position {fields[0]}:{fields[1]}. Expected phased data with '|' separator.\n")
                        continue
                    
                    # Clear INFO field if requested
                    if not keep_info:
                        fields[7] = "."
                    
                    # Process FORMAT fields
                    if keep_format:
                        # Keep all FORMAT fields, but update GT
                        gt_fields[0] = new_gt
                        new_gt_data = ':'.join(gt_fields)
                    else:
                        # Only keep GT field
                        new_gt_data = new_gt
                        fields[8] = "GT"  # Update FORMAT to only include GT
                    
                    # Create new line with just the first 9 columns + modified sample
                    new_line = fields[:9] + [new_gt_data]
                    output_line = '\t'.join(new_line) + '\n'
                    
                    # Write in the appropriate format
                    if output_is_gzipped:
                        out_file.write(output_line.encode('utf-8'))
                    else:
                        out_file.write(output_line)
        finally:
            out_file.close()

def main():
    args = parse_arguments()
    
    # Check if input file exists
    if not os.path.exists(args.input_vcf):
        sys.stderr.write(f"Error: Input file '{args.input_vcf}' does not exist.\n")
        sys.exit(1)
    
    # Process the VCF file
    try:
        process_vcf(args.input_vcf, args.output_vcf, args.sample, args.haplotype, 
                   args.keep_format, args.keep_info)
        print(f"Extracted haplotype {args.haplotype} from sample {args.sample}")
        print(f"Output written to {args.output_vcf}")
    except TypeError as e:
        print(f"Error: {e}")
        print("This might be due to mixing binary and text modes in file operations.")
        sys.exit(1)

if __name__ == "__main__":
    main()
