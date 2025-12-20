#!/usr/bin/env python3
"""
Batch generate iTOL-compatible annotation files for updated genus statistics.

This script processes genus_stats_merged_*.csv files with new columns and generates:
1. Bar chart for Genus_Number
2. Bar chart for Mean_Plasmid_Percent
3. Heatmap for detailed Defense types (Per_kb)
4. Heatmap for detailed AntiDS types (Per_kb)
5. Heatmap for detailed AMR types (Per_kb)
6. Heatmap for detailed Antidefense types (Per_kb) - specific antidefense columns
7. Class-based color annotation

Usage:
    python 206_batch_generate_itol_annotations.py -i <input_dir> -o <output_dir> [-h]

Arguments:
    -i: Input directory containing genus_stats_merged_*.csv files
    -o: Output directory for annotation files
    -h: Show this help message

Example:
    python 206_batch_generate_itol_annotations.py -i Result/NCBI_4395_Batch/05_Tree/Genus_Level/Tree_annotation_file_prepare -o Result/NCBI_4395_Batch/05_Tree/Genus_Level/Tree_visualization
"""

import pandas as pd
import argparse
import numpy as np
import os
import glob
from pathlib import Path

# Class color mapping
CLASS_COLORS = {
    'Actinomycetia': '#98df8a',
    'Alphaproteobacteria': '#aec7e8',
    'Bacilli': '#ff7f0e',
    'Thermoleophilia': '#ff9896',
    'Bacteroidia': '#d62728',
    'Gammaproteobacteria': '#ffbb78',
    'Deinococci': '#1f77b4',
    'Acidimicrobiia': '#2ca02c',
    'Campylobacteria': '#9467bd'
}


def generate_bar_chart(df, column, output_file, label, color='#5499C7'):
    """
    Generate iTOL bar chart annotation file.
    
    Args:
        df: Input DataFrame
        column: Column name for bar values
        output_file: Output file path
        label: Dataset label
        color: Bar color (hex)
    """
    lines = [
        "DATASET_SIMPLEBAR",
        "SEPARATOR COMMA",
        f"DATASET_LABEL,{label}",
        f"COLOR,{color}",
        "",
        f"FIELD_LABELS,{label}",
        f"FIELD_COLORS,{color}",
        f"LEGEND_TITLE,{label}",
        "LEGEND_SHAPES,1",
        f"LEGEND_COLORS,{color}",
        f"LEGEND_LABELS,{label}",
        "",
        "DATA"
    ]
    
    for _, row in df.iterrows():
        genus = row['Genus']
        value = row[column]
        lines.append(f"{genus},{value}")
    
    with open(output_file, "w") as f:
        f.write("\n".join(lines))


def generate_multi_column_heatmap(df, columns, output_file, label, colors=None, use_log=True):
    """
    Generate iTOL heatmap annotation file for multiple columns.
    
    Args:
        df: Input DataFrame
        columns: List of column names for heatmap
        output_file: Output file path
        label: Dataset label
        colors: List of colors for each column (hex), or single color
        use_log: Whether to apply log10 transformation
    """
    # If colors is a single string, use it for all columns
    if isinstance(colors, str):
        colors = [colors] * len(columns)
    elif colors is None:
        colors = ['#E74C3C'] * len(columns)
    
    # Clean column names for labels (remove _Per_kb suffix)
    field_labels = [col.replace('_Per_kb', '') for col in columns]
    
    legend_title = f"log10(Per kb + 1)" if use_log else "Per kb"
    
    lines = [
        "DATASET_HEATMAP",
        "SEPARATOR COMMA",
        f"DATASET_LABEL,{label}",
        f"COLOR,{colors[0]}",
        "",
        f"FIELD_LABELS,{','.join(field_labels)}",
        f"FIELD_COLORS,{','.join(colors)}",
        f"LEGEND_TITLE,{legend_title}",
        "LEGEND_SHAPES," + ",".join(['1'] * len(columns)),
        f"LEGEND_COLORS,{','.join(colors)}",
        f"LEGEND_LABELS,{','.join(field_labels)}",
        "",
        "DATA"
    ]
    
    for _, row in df.iterrows():
        genus = row['Genus']
        values = []
        for col in columns:
            value = row[col] if col in row else 0
            # Apply log10 transformation if requested
            if use_log:
                transformed_value = np.log10(value + 1)
            else:
                transformed_value = value
            values.append(str(transformed_value))
        lines.append(f"{genus},{','.join(values)}")
    
    with open(output_file, "w") as f:
        f.write("\n".join(lines))


def generate_class_colors(df, output_file):
    """
    Generate iTOL colorstrip annotation file based on Class column.
    
    Args:
        df: Input DataFrame with 'Genus' and 'Class' columns
        output_file: Output file path
    """
    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR COMMA",
        "DATASET_LABEL,Class",
        "COLOR,#000000",
        "",
        "LEGEND_TITLE,Class",
        "LEGEND_SHAPES," + ",".join(['1'] * len(CLASS_COLORS)),
        "LEGEND_COLORS," + ",".join(CLASS_COLORS.values()),
        "LEGEND_LABELS," + ",".join(CLASS_COLORS.keys()),
        "",
        "DATA"
    ]
    
    for _, row in df.iterrows():
        genus = row['Genus']
        class_name = row['Class'] if 'Class' in row and pd.notna(row['Class']) else 'Unknown'
        color = CLASS_COLORS.get(class_name, '#cccccc')  # Default gray for unknown
        lines.append(f"{genus},{color},{class_name}")
    
    with open(output_file, "w") as f:
        f.write("\n".join(lines))


def process_file(input_file, output_dir):
    """
    Process a single genus_stats_merged file and generate all annotation files.
    
    Args:
        input_file: Path to input CSV file
        output_dir: Directory for output files
        
    Returns:
        Tuple of (success, message, annotation_count)
    """
    # Extract prefix from filename
    basename = os.path.basename(input_file)
    # Remove 'genus_stats_merged_' prefix and '.csv' suffix
    prefix = basename.replace('genus_stats_merged_', '').replace('.csv', '')
    
    try:
        # Load data
        df = pd.read_csv(input_file)
        
        # Filter out specific genera
        exclude_genera = ['Aliarcobacter', 'Saccharopolyspora', 'Ralstonia', 'Blastococcus']
        df = df[~df['Genus'].isin(exclude_genera)]
        
        # Check if required columns exist
        required_cols = ['Genus', 'Genus_Number', 'Mean_Plasmid_Percent', 'Class']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            return False, f"Missing columns: {missing_cols}", 0
        
        annotation_count = 0
        
        # 1. Generate bar chart for Genus_Number
        bar_genus_output = os.path.join(output_dir, f"{prefix}_genus_number_bar.txt")
        generate_bar_chart(
            df, 
            'Genus_Number', 
            bar_genus_output, 
            'Genus Number',
            color='#3498DB'  # Blue
        )
        annotation_count += 1
        
        # 2. Generate bar chart for Mean_Plasmid_Percent
        bar_plasmid_output = os.path.join(output_dir, f"{prefix}_plasmid_percent_bar.txt")
        generate_bar_chart(
            df, 
            'Mean_Plasmid_Percent', 
            bar_plasmid_output, 
            'Mean Plasmid Percent',
            color='#2ECC71'  # Green
        )
        annotation_count += 1
        
        # 3. Generate Class color strip
        class_output = os.path.join(output_dir, f"{prefix}_class_colors.txt")
        generate_class_colors(df, class_output)
        annotation_count += 1
        
        # 4. Generate heatmap for Defense types (Per_kb)
        defense_cols = ['PD-T4-6_Per_kb', 'RM type II_Per_kb', 'SoFic_Per_kb', 
                       'RM type I_Per_kb', 'RM type IV_Per_kb', 'Ceres_Per_kb',
                       'AbiE_Per_kb', 'Wadjet type I_Per_kb', 'Gabija_Per_kb',
                       'CRISPR-Cas type I_Per_kb', 'Zorya type III_Per_kb', 
                       'Septu_Per_kb', 'Defense_Others_Per_kb']
        # Filter columns that exist in dataframe
        defense_cols_exist = [col for col in defense_cols if col in df.columns]
        
        if defense_cols_exist:
            defense_output = os.path.join(output_dir, f"{prefix}_defense_per_kb_heatmap.txt")
            generate_multi_column_heatmap(
                df,
                defense_cols_exist,
                defense_output,
                'Defense Systems (Per kb)',
                colors='#E74C3C',  # Red
                use_log=True
            )
            annotation_count += 1
        
        # 5. Generate heatmap for AntiDS types (Per_kb)
        antids_cols = [col for col in df.columns if col.startswith('Anti_') and col.endswith('_Per_kb')]
        
        if antids_cols:
            antids_output = os.path.join(output_dir, f"{prefix}_antids_per_kb_heatmap.txt")
            generate_multi_column_heatmap(
                df,
                antids_cols,
                antids_output,
                'Anti-Defense Systems (Per kb)',
                colors='#9B59B6',  # Purple
                use_log=True
            )
            annotation_count += 1
        
        # 6. Generate heatmap for detailed Antidefense types (Per_kb)
        # Specific antidefense columns including NADP and Other
        antidefense_detail_cols = [
            'Anti_CBASS_Per_kb', 'Anti_CRISPR_Per_kb', 'Anti_Dnd_Per_kb',
            'Anti_Gabija_Per_kb', 'Anti_Pycsar_Per_kb', 'Anti_RM_Per_kb',
            'Anti_Thoeris_Per_kb', 'NADP_Per_kb', 'Other_Per_kb'
        ]
        # Filter columns that exist in dataframe
        antidefense_detail_cols_exist = [col for col in antidefense_detail_cols if col in df.columns]
        
        if antidefense_detail_cols_exist:
            antidefense_detail_output = os.path.join(output_dir, f"{prefix}_antidefense_detail_per_kb_heatmap.txt")
            generate_multi_column_heatmap(
                df,
                antidefense_detail_cols_exist,
                antidefense_detail_output,
                'Antidefense Detail (Per kb)',
                colors='#8E44AD',  # Dark purple
                use_log=True
            )
            annotation_count += 1
        
        # 7. Generate heatmap for AMR types (Per_kb)
        # Get only AMR-related columns (excluding Defense and Anti-Defense)
        amr_per_kb_cols = []
        for col in df.columns:
            if col.endswith('_Per_kb'):
                col_base = col.replace('_Per_kb', '')
                # Check if it's an AMR column (contains common AMR gene names)
                # Exclude Defense-related and Anti-Defense columns
                if (any(amr_name in col_base for amr_name in ['bla', 'arr', 'ampC', 'cml', 'aac', 'erm', 
                                                               'aph', 'vanR', 'catB', 'fosX', 'rox', 
                                                               'cpt', 'vat', 'vanX', 'catA', 'tet', 
                                                               'lsa', 'AMR_Others'])
                    and not col.startswith('Anti_')
                    and 'Defense' not in col
                    and col not in defense_cols_exist):
                    amr_per_kb_cols.append(col)
        
        if amr_per_kb_cols:
            amr_output = os.path.join(output_dir, f"{prefix}_amr_per_kb_heatmap.txt")
            generate_multi_column_heatmap(
                df,
                amr_per_kb_cols,
                amr_output,
                'AMR Genes (Per kb)',
                colors='#F39C12',  # Orange
                use_log=True
            )
            annotation_count += 1
        
        return True, f"Successfully processed {prefix} ({len(df)} genera, {annotation_count} annotations)", annotation_count
        
    except Exception as e:
        import traceback
        return False, f"Error processing {prefix}: {str(e)}\n{traceback.format_exc()}", 0


def main():
    """Main function to batch process all genus_stats_merged files."""
    parser = argparse.ArgumentParser(
        description="Batch generate iTOL annotation files (updated version)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '-i', '--input_dir',
        required=True,
        help='Input directory containing genus_stats_merged_*.csv files'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=True,
        help='Output directory for annotation files'
    )
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Find all genus_stats_merged_*.csv files
    pattern = os.path.join(args.input_dir, 'genus_stats_merged_*.csv')
    input_files = sorted(glob.glob(pattern))
    
    if not input_files:
        print(f"Error: No genus_stats_merged_*.csv files found in {args.input_dir}")
        return
    
    print("=" * 80)
    print("Batch Processing iTOL Annotation Files (Updated Version)")
    print("=" * 80)
    print(f"\nInput directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Found {len(input_files)} files to process\n")
    
    # Process each file
    success_count = 0
    fail_count = 0
    total_annotations = 0
    
    for i, input_file in enumerate(input_files, 1):
        print(f"[{i}/{len(input_files)}] Processing {os.path.basename(input_file)}...")
        success, message, annotation_count = process_file(input_file, args.output_dir)
        
        if success:
            print(f"  ✓ {message}")
            success_count += 1
            total_annotations += annotation_count
        else:
            print(f"  ✗ {message}")
            fail_count += 1
    
    print("\n" + "=" * 80)
    print("Batch Processing Complete")
    print("=" * 80)
    print(f"\nTotal files processed: {len(input_files)}")
    print(f"  Successful: {success_count}")
    print(f"  Failed: {fail_count}")
    print(f"\nTotal annotation files generated: {total_annotations}")
    print(f"  - Genus Number bar charts: {success_count}")
    print(f"  - Plasmid Percent bar charts: {success_count}")
    print(f"  - Class color strips: {success_count}")
    print(f"  - Defense Per kb heatmaps: {success_count}")
    print(f"  - AntiDS Per kb heatmaps: {success_count}")
    print(f"  - Antidefense Detail Per kb heatmaps: {success_count}")
    print(f"  - AMR Per kb heatmaps: {success_count}")
    print("\nOutput location: " + args.output_dir)
    print("=" * 80)


if __name__ == "__main__":
    main()

