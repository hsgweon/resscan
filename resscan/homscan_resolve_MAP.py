#!/usr/bin/env python3

# resscan/homscan_resolve_MAP.py
import argparse
import sys
import re
import csv
from collections import defaultdict
import pandas as pd
import numpy as np

# Formatting settings
pd.set_option('display.width', 120)
pd.set_option('display.max_columns', 15)
pd.set_option('display.max_colwidth', 40)

class BColors:
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
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
    
    if metric_col in df.columns:
        df[metric_col] = pd.to_numeric(df[metric_col], errors='coerce').fillna(0)
    else:
        print(BColors.red(f"Error: Metric column '{metric_col}' not found."), file=sys.stderr)
        sys.exit(1)

    all_aros = set()
    aro_pattern = re.compile(r'ARO_\d+')
    
    # Check the new 'ARO' column for identifiers
    for aro_val in df['ARO'].astype(str):
        aros_in_key = aro_pattern.findall(aro_val)
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
        print(BColors.yellow(f"Warning: Priors file not found at '{filepath}'."), file=sys.stderr)
        return {}
    return priors

def run_iterative_solver(detailed_df, all_aros_list, priors, base_prior, prior_strength, metric_col: str, resolved_metric_col: str):
    print(BColors.cyan(f"--- Running solver using metric: {metric_col} ---"))
    
    calc_col = metric_col
    if detailed_df[metric_col].sum() == 0:
        for fallback in ['RPKG', 'RPK', 'Read_Count']:
            if fallback in detailed_df.columns and detailed_df[fallback].sum() > 0:
                calc_col = fallback
                break
    
    aro_pattern = re.compile(r'ARO_\d+')
    unique_abundances = defaultdict(float)
    
    # Separate ambiguous vs unique based on the 'ARO' column
    for _, row in detailed_df.iterrows():
        aro_val = str(row['ARO'])
        if 'multiple' not in aro_val:
            aros = aro_pattern.findall(aro_val)
            if aros:
                unique_abundances[aros[0]] += row[calc_col]
    
    static_abundances = pd.Series(unique_abundances, dtype=float).reindex(all_aros_list, fill_value=0.0)
    for aro, hist_count in priors.items():
        if aro in static_abundances.index:
            static_abundances[aro] += hist_count * prior_strength * 0.1 
    
    static_abundances += base_prior * 1e-6
    ambiguous_rows = detailed_df[detailed_df['ARO'].str.contains('multiple', na=False)].copy()
    
    # If a row is multiple, we must find all AROs associated with that family
    # Here we look at the ARO column content which now looks like multiple[ARO_1|ARO_2]
    ambiguous_rows['AROs_in_group'] = ambiguous_rows['ARO'].apply(
        lambda x: sorted(list(set(aro_pattern.findall(x))))
    )

    final_abundances = static_abundances.copy()
    for i in range(50):
        abundances_before_iter = final_abundances.copy()
        ambiguous_abun_this_iter = pd.Series(0.0, index=all_aros_list)
        for _, row in ambiguous_rows.iterrows():
            aros_in_group = row['AROs_in_group']
            if not aros_in_group: continue
            current_abun = abundances_before_iter.loc[aros_in_group]
            total_w = current_abun.sum()
            props = current_abun / total_w if total_w > 1e-9 else pd.Series(1.0/len(aros_in_group), index=aros_in_group)
            ambiguous_abun_this_iter = ambiguous_abun_this_iter.add(row[calc_col] * props, fill_value=0)
        
        final_abundances = static_abundances.add(ambiguous_abun_this_iter, fill_value=0)
        if np.allclose(final_abundances.values, abundances_before_iter.values, atol=1e-6, rtol=1e-5):
            break

    final_resolved_abundances = pd.Series(0.0, index=all_aros_list)
    for _, row in detailed_df.iterrows():
        aro_val = str(row['ARO'])
        if 'multiple' not in aro_val:
            aros = aro_pattern.findall(aro_val)
            if aros: final_resolved_abundances[aros[0]] += row[metric_col]

    for _, row in ambiguous_rows.iterrows():
        aros_in_group = row['AROs_in_group']
        if not aros_in_group: continue
        curr_w = final_abundances.loc[aros_in_group]
        props = curr_w / curr_w.sum() if curr_w.sum() > 1e-9 else pd.Series(1.0/len(aros_in_group), index=aros_in_group)
        final_resolved_abundances = final_resolved_abundances.add(row[metric_col] * props, fill_value=0)

    return final_resolved_abundances.reset_index().rename(columns={'index':'ARO', 0:resolved_metric_col}), \
           final_abundances.reset_index().rename(columns={'index':'ARO', 0:'Weight'})

def generate_final_summary(resolved_df, weight_df, detailed_df, metric_col, resolved_metric_col):
    print(BColors.cyan("--- Generating final summary report... ---"))
    aro_pattern = re.compile(r'ARO_\d+')
    weights_map = dict(zip(weight_df['ARO'], weight_df['Weight']))
    
    summary_rows = []
    processed_families = set()

    for _, row in detailed_df.iterrows():
        family = row['AMR_Gene_Family']
        aro_val = str(row['ARO'])

        if 'multiple' not in aro_val:
            new_row = row.to_dict()
            aros = aro_pattern.findall(aro_val)
            if aros:
                aro = aros[0]
                new_row['Top_ARO'] = aro
                new_row['Allocation_Proportions'] = f"{aro}:100.0%"
                new_row[resolved_metric_col] = resolved_df.loc[resolved_df['ARO'] == aro, resolved_metric_col].values[0]
                summary_rows.append(new_row)
        
        elif family not in processed_families:
            family_rows = detailed_df[detailed_df['AMR_Gene_Family'] == family]
            all_aros_for_fam = sorted(list(set(aro_pattern.findall(" ".join(family_rows['ARO'].astype(str))))))
            
            agg_row = {col: 'multiple' for col in detailed_df.columns}
            agg_row['AMR_Gene_Family'] = family
            agg_row['ARO'] = "multiple"
            
            sum_cols = ['Read_Count', 'Fragment_Count']
            num_cols = family_rows.select_dtypes(include=np.number).columns
            agg_row.update(family_rows[sum_cols].sum().to_dict())
            agg_row.update(family_rows[[c for c in num_cols if c not in sum_cols]].mean().to_dict())

            fam_weights = {a: weights_map.get(a, 0.0) for a in all_aros_for_fam}
            total_w = sum(fam_weights.values())
            agg_row[resolved_metric_col] = resolved_df[resolved_df['ARO'].isin(all_aros_for_fam)][resolved_metric_col].sum()

            if total_w > 1e-9:
                props = {a: w/total_w for a, w in fam_weights.items() if w > 0}
                winner = max(props, key=props.get)
                agg_row['Top_ARO'] = winner if list(props.values()).count(props[winner]) == 1 else 'N/A'
                sorted_p = sorted(props.items(), key=lambda x: x[1], reverse=True)
                agg_row['Allocation_Proportions'] = ",".join([f"{a}:{p*100:.1f}%" for a, p in sorted_p])
            else:
                agg_row['Top_ARO'], agg_row['Allocation_Proportions'] = "None", "No Evidence"
            
            summary_rows.append(agg_row)
            processed_families.add(family)

    return pd.DataFrame(summary_rows)

def main():
    parser = argparse.ArgumentParser(description="Resolves ambiguous gene families via MAP.")
    parser.add_argument("-i", "--input-file", required=True)
    parser.add_argument("-p", "--priors-file")
    parser.add_argument("-o", "--output-prefix", required=True)
    parser.add_argument("--metric-column", default="RPKG")
    parser.add_argument("--base-prior", type=float, default=1.0)
    parser.add_argument("--prior-strength", type=float, default=1.0)
    args = parser.parse_args()

    output_filename = f"{args.output_prefix}_homscan_MAP.tsv"
    resolved_metric_col = f"Resolved_{args.metric_column}"

    detailed_df, all_aros_list = parse_input_table(args.input_file, args.metric_column)
    priors = parse_priors_file(args.priors_file)
    
    if detailed_df.empty:
        sys.exit(0)

    resolved_df, weight_df = run_iterative_solver(detailed_df, all_aros_list, priors, args.base_prior, args.prior_strength, args.metric_column, resolved_metric_col)
    summary_df = generate_final_summary(resolved_df, weight_df, detailed_df, args.metric_column, resolved_metric_col)
    
    cols_to_keep = [
        'AMR_Gene_Family', 'ARO', 'Read_Count', 'Fragment_Count', 
        'Lateral_Coverage_%', 'Gene_Length_bp', 'RPK', 'FPK', 'RPKG', 'FPKG', 
        'RPKPC', 'FPKPC', 'RPKPMC', 'FPKPMC', resolved_metric_col, 'Top_ARO', 'Allocation_Proportions'
    ]
    
    for col in cols_to_keep:
        if col not in summary_df.columns: summary_df[col] = "N/A"

    summary_df.reindex(columns=cols_to_keep).to_csv(output_filename, sep='\t', index=False, float_format='%.6f')
    print(BColors.green(f"--- Successfully wrote final summary to: {output_filename} ---"))

if __name__ == "__main__":
    main()