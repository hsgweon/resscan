#!/usr/bin/env python3

# resscan/homscan_resolve_MAP.py
import argparse
import sys
import re
import json
import csv
from collections import defaultdict

import pandas as pd
import numpy as np

pd.set_option('display.width', 120)
pd.set_option('display.max_columns', 15)
pd.set_option('display.max_colwidth', 40)

class BColors:
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    enabled = sys.stdout.isatty()
    @classmethod
    def _colorize(cls, c, t): return f"{c}{t}{cls.ENDC}" if cls.enabled else t
    @classmethod
    def cyan(cls, t): return cls._colorize(cls.OKCYAN, t)
    @classmethod
    def green(cls, t): return cls._colorize(cls.OKGREEN, t)
    @classmethod
    def red(cls, t): return cls._colorize(cls.FAIL, t)
    @classmethod
    def yellow(cls, t): return cls._colorize(cls.WARNING, t)

def parse_input_table(filepath: str, metric_col: str):
    print(BColors.cyan(f"--- Parsing input data from: {filepath} ---"))
    try:
        df = pd.read_csv(filepath, sep='\t')
    except FileNotFoundError:
        print(BColors.red(f"Error: Input file not found at '{filepath}'"), file=sys.stderr)
        sys.exit(1)
    
    if metric_col not in df.columns:
        # Check if the column is present but 0 due to N/A filling in previous steps
        pass 

    # Ensure metric column is numeric
    if metric_col in df.columns:
        df[metric_col] = pd.to_numeric(df[metric_col], errors='coerce').fillna(0)
    else:
        print(BColors.red(f"Error: Metric column '{metric_col}' not found."), file=sys.stderr)
        sys.exit(1)

    all_aros = set()
    aro_pattern = re.compile(r'ARO_\d+')
    
    for key in df['#AMR_Gene_Family;ARO']:
        aros_in_key = aro_pattern.findall(key)
        all_aros.update(aros_in_key)

    all_aros_list = sorted(list(all_aros))
    
    print(BColors.green(f"--- Found {len(df)} total observations for {len(all_aros_list)} unique AROs."))
    return df, all_aros_list

def parse_priors_file(filepath: str):
    if not filepath: return {}
    print(BColors.cyan(f"--- Loading priors from: {filepath} ---"))
    priors = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.strip().split('\t')
                if len(parts) >= 1:
                    priors[parts[0]] = float(parts[1]) if len(parts) > 1 else 1.0
    except FileNotFoundError:
        print(BColors.yellow(f"Warning: Priors file not found at '{filepath}'. Proceeding without explicit priors."), file=sys.stderr)
        return {}
    print(BColors.green(f"--- Loaded {len(priors)} historical/clinical priors."))
    return priors

def run_iterative_solver(detailed_df, all_aros_list, priors, base_prior, prior_strength, metric_col: str, resolved_metric_col: str):
    """
    Resolves ambiguous abundances using a fast, iterative allocation algorithm (EM-like).
    This finds the Maximum A Posteriori (MAP) estimate.
    """
    print(BColors.cyan(f"--- Running solver using metric: {metric_col} ---"))
    
    # 1. Check if the metric has data. If not, fallback for PROPORTION calculation only.
    calc_col = metric_col
    total_abun = detailed_df[metric_col].sum()
    
    if total_abun == 0:
        print(BColors.yellow(f"Warning: Selected metric '{metric_col}' sums to 0 (likely due to missing SCGs)."))
        # Try fallbacks for calculation logic
        for fallback in ['RPKG', 'RPK', 'Read_Count']:
            if fallback in detailed_df.columns and detailed_df[fallback].sum() > 0:
                print(BColors.yellow(f"-> Falling back to '{fallback}' to determine resolution probabilities."))
                calc_col = fallback
                break
    
    aro_pattern = re.compile(r'ARO_\d+')

    # 2. Initial Estimation (Unique Reads)
    unique_abundances = defaultdict(float)
    for _, row in detailed_df.iterrows():
        key = row['#AMR_Gene_Family;ARO']
        if 'multiple[' not in key and not key.endswith(';multiple'):
            aros = aro_pattern.findall(key)
            if aros:
                unique_abundances[aros[0]] += row[calc_col]
    
    static_abundances = pd.Series(unique_abundances, dtype=float).reindex(all_aros_list, fill_value=0.0)

    # 3. Apply Priors
    for aro, hist_count in priors.items():
        if aro in static_abundances.index:
            static_abundances[aro] += hist_count * prior_strength * 0.1 
    
    static_abundances += base_prior * 1e-6

    ambiguous_rows = detailed_df[
        detailed_df['#AMR_Gene_Family;ARO'].str.contains('multiple', na=False)
    ].copy()
    
    ambiguous_rows['AROs_in_group'] = ambiguous_rows['#AMR_Gene_Family;ARO'].apply(
        lambda x: sorted(list(set(aro_pattern.findall(x))))
    )

    final_abundances = static_abundances.copy()

    print(BColors.cyan("--- Starting iterative allocation... ---"))
    for i in range(50):
        abundances_before_iter = final_abundances.copy()
        ambiguous_abun_this_iter = pd.Series(0.0, index=all_aros_list)

        for _, row in ambiguous_rows.iterrows():
            aros_in_group = row['AROs_in_group']
            if not aros_in_group: continue

            current_abun_of_group_members = abundances_before_iter.loc[aros_in_group]
            total_group_weight = current_abun_of_group_members.sum()

            if total_group_weight > 1e-9:
                proportions = current_abun_of_group_members / total_group_weight
            else:
                num_members = len(aros_in_group)
                if num_members > 0:
                    proportions = pd.Series(1.0 / num_members, index=aros_in_group)
                else:
                    continue

            # Distribute the *Calculation* Metric to update weights for next round
            distributed_calc_abun = row[calc_col] * proportions
            ambiguous_abun_this_iter = ambiguous_abun_this_iter.add(distributed_calc_abun, fill_value=0)

        final_abundances = static_abundances.add(ambiguous_abun_this_iter, fill_value=0)
        
        if np.allclose(final_abundances.values, abundances_before_iter.values, atol=1e-6, rtol=1e-5):
            print(BColors.green(f"--- Converged after {i+1} iterations. ---"))
            break
    else:
        print(BColors.yellow("--- Reached maximum iterations without convergence. ---"))

    # 4. Final Pass: Apply calculated proportions to the ORIGINAL requested metric
    # This ensures the output column 'Resolved_FPKPMC' matches the input 'FPKPMC' (even if it's 0)
    final_resolved_abundances = pd.Series(0.0, index=all_aros_list)
    
    # Add unique
    for _, row in detailed_df.iterrows():
        key = row['#AMR_Gene_Family;ARO']
        if 'multiple[' not in key and not key.endswith(';multiple'):
            aros = aro_pattern.findall(key)
            if aros:
                final_resolved_abundances[aros[0]] += row[metric_col]

    # Add ambiguous using final weights
    for _, row in ambiguous_rows.iterrows():
        aros_in_group = row['AROs_in_group']
        if not aros_in_group: continue
        
        current_weights = final_abundances.loc[aros_in_group]
        total_weight = current_weights.sum()
        
        if total_weight > 1e-9:
            props = current_weights / total_weight
        else:
            props = pd.Series(1.0 / len(aros_in_group), index=aros_in_group)
            
        distributed = row[metric_col] * props
        final_resolved_abundances = final_resolved_abundances.add(distributed, fill_value=0)

    # Return DataFrame for merging
    result_df = final_resolved_abundances.reset_index()
    result_df.columns = ['ARO', resolved_metric_col]
    
    # Also return the final_abundances (weights) so we can calculate Top_ARO correctly 
    # even if resolved_metric is 0
    weight_df = final_abundances.reset_index()
    weight_df.columns = ['ARO', 'Weight']
    
    return result_df, weight_df

def generate_final_summary(resolved_df, weight_df, detailed_df, metric_col: str, resolved_metric_col: str):
    print(BColors.cyan("--- Generating final summary report... ---"))
    
    aro_pattern = re.compile(r'ARO_\d+')
    
    # Map weights for quick lookup
    weights_map = dict(zip(weight_df['ARO'], weight_df['Weight']))

    family_to_ambiguous_aros = defaultdict(set)
    for _, row in detailed_df.iterrows():
        if 'multiple[' in row['#AMR_Gene_Family;ARO'] or row['#AMR_Gene_Family;ARO'].endswith(';multiple'):
            family = row['#AMR_Gene_Family;ARO'].split(';')[0]
            aros = aro_pattern.findall(row['#AMR_Gene_Family;ARO'])
            family_to_ambiguous_aros[family].update(aros)

    summary_rows = []
    processed_families = set()

    for index, row in detailed_df.iterrows():
        key = row['#AMR_Gene_Family;ARO']
        family = key.split(';')[0]

        if 'multiple[' not in key and not key.endswith(';multiple'):
            new_row = row.to_dict()
            aros = aro_pattern.findall(key)
            if aros:
                aro = aros[0]
                new_row['Top_ARO'] = aro
                new_row['Allocation_Proportions'] = f"{aro}:100.0%"
                # Add resolved metric
                new_row[resolved_metric_col] = resolved_df.loc[resolved_df['ARO'] == aro, resolved_metric_col].values[0]
                summary_rows.append(new_row)

        elif family not in processed_families:
            ambiguous_aros_for_family = list(family_to_ambiguous_aros[family])
            
            agg_row = {col: 'multiple' for col in detailed_df.columns}
            agg_row['#AMR_Gene_Family;ARO'] = f"{family};multiple"
            
            family_multiple_rows = detailed_df[
                (detailed_df['#AMR_Gene_Family;ARO'].str.startswith(f"{family};multiple"))
            ]
            
            cols_to_sum = ['Read_Count', 'Fragment_Count']
            all_numeric_cols = family_multiple_rows.select_dtypes(include=np.number).columns
            cols_to_average = [col for col in all_numeric_cols if col not in cols_to_sum]

            summed_vals = family_multiple_rows[cols_to_sum].sum().to_dict()
            averaged_vals = family_multiple_rows[cols_to_average].mean().to_dict()
            
            agg_row.update(summed_vals)
            agg_row.update(averaged_vals)

            # Determine proportions based on WEIGHTS (calculated via fallback), not necessarily the resolved metric (which might be 0)
            family_weights = {aro: weights_map.get(aro, 0.0) for aro in ambiguous_aros_for_family}
            total_weight = sum(family_weights.values())
            
            # Calculate Resolved Metric Sum for the family
            family_resolved_metric_sum = resolved_df[resolved_df['ARO'].isin(ambiguous_aros_for_family)][resolved_metric_col].sum()
            agg_row[resolved_metric_col] = family_resolved_metric_sum

            if total_weight > 1e-9:
                proportions = {aro: w / total_weight for aro, w in family_weights.items() if w > 0}
                
                # Determine Top ARO
                winner_aro = max(proportions, key=proportions.get)
                # Check for tie
                max_val = proportions[winner_aro]
                if list(proportions.values()).count(max_val) > 1:
                     agg_row['Top_ARO'] = 'N/A'
                else:
                     agg_row['Top_ARO'] = winner_aro
                
                sorted_proportions = sorted(proportions.items(), key=lambda item: item[1], reverse=True)
                proportions_str = ",".join([f"{aro}:{prop*100:.1f}%" for aro, prop in sorted_proportions])
                agg_row['Allocation_Proportions'] = proportions_str
            else:
                agg_row['Top_ARO'] = "None"
                agg_row['Allocation_Proportions'] = "No Evidence"
            
            summary_rows.append(agg_row)
            processed_families.add(family)

    summary_df = pd.DataFrame(summary_rows)
    return summary_df

def main():
    parser = argparse.ArgumentParser(
        description="Resolves ambiguous gene families via MAP."
    )
    parser.add_argument("-i", "--input-file", required=True, help="Input detailed report.")
    parser.add_argument("-p", "--priors-file", help="Priors file.")
    parser.add_argument("-o", "--output-prefix", required=True, help="Output prefix.")
    parser.add_argument("--metric-column", default="RPKG", help="Metric column.")
    parser.add_argument("--base-prior", type=float, default=1.0)
    parser.add_argument("--prior-strength", type=float, default=1.0)
    args = parser.parse_args()

    output_filename = f"{args.output_prefix}_homscan_MAP.tsv"
    resolved_metric_col = f"Resolved_{args.metric_column}"

    detailed_df, all_aros_list = parse_input_table(args.input_file, args.metric_column)
    priors = parse_priors_file(args.priors_file)
    
    if detailed_df.empty:
        print(BColors.yellow("Warning: Input file is empty."), file=sys.stderr)
        with open(output_filename, 'w') as f: f.write("")
        sys.exit(0)

    # Run solver
    resolved_df, weight_df = run_iterative_solver(detailed_df, all_aros_list, priors, args.base_prior, args.prior_strength, args.metric_column, resolved_metric_col)
    
    # Generate Summary
    summary_df = generate_final_summary(resolved_df, weight_df, detailed_df, args.metric_column, resolved_metric_col)
    
    try:
        # Organize columns
        cols_to_keep = [
            '#AMR_Gene_Family;ARO', 'Read_Count', 'Fragment_Count', 
            'Lateral_Coverage_%', 'Gene_Length_bp', 'RPK', 'FPK', 'RPKG', 'FPKG', 
            'RPKPC', 'FPKPC', 'RPKPMC', 'FPKPMC', resolved_metric_col, 'Top_ARO', 'Allocation_Proportions'
        ]
        
        # Ensure cols exist
        for col in cols_to_keep:
            if col not in summary_df.columns: summary_df[col] = "N/A"

        summary_df = summary_df.reindex(columns=cols_to_keep)
        
        summary_df.to_csv(output_filename, sep='\t', index=False, float_format='%.6f')
        print(BColors.green(f"\n--- Successfully wrote final summary to: {output_filename} ---"))

    except IOError as e:
        print(BColors.red(f"Error writing output file: {e}"), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()