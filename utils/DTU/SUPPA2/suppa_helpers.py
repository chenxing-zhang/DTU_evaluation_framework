# suppa_helpers.py
import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import mannwhitneyu

def psiPerIsoform_direct(tpm_df, ioi_df):
    """
    Calculates PSI values for each transcript directly from TPM and IOI dataframes.

    Args:
        tpm_df (pd.DataFrame): Dataframe with transcripts as index and samples as columns, containing TPM values.
        ioi_df (pd.DataFrame): Dataframe with at least 'transcript' and 'gene' columns.

    Returns:
        pd.DataFrame: Dataframe with transcripts as index and samples as columns, containing PSI values.
                      Returns None if input is invalid or calculation fails.
    """
    try:
        # Ensure ioi_df has required columns
        if not all(col in ioi_df.columns for col in ['transcript', 'gene']):
            print("Error: ioi_df must contain 'transcript' and 'gene' columns.")
            return None

        # Ensure tpm_df index is named or reset it to use as a column
        if tpm_df.index.name is None:
            tpm_df.index.name = 'transcript'

        # Make sure index is string type if it's not already, for merging
        tpm_df.index = tpm_df.index.astype(str)
        ioi_df['transcript'] = ioi_df['transcript'].astype(str)

        # Reset index to merge tpm data with gene info
        tpm_reset = tpm_df.reset_index()

        # Merge TPM data with IOI data to get gene information
        # Keep only transcripts present in both tpm_df and ioi_df
        merged_df = pd.merge(tpm_reset, ioi_df[['transcript', 'gene']], on='transcript', how='inner')

        if merged_df.empty:
            print("Error: No common transcripts found between TPM data and IOI data.")
            return None

        # Set index back to transcript for easier processing (optional, could work with multi-index)
        # We need gene and transcript info for grouping
        # merged_df = merged_df.set_index('transcript') # Keep transcript column for grouping

        # Identify sample columns (all columns except 'transcript' and 'gene')
        sample_cols = [col for col in merged_df.columns if col not in ['transcript', 'gene']]

        if not sample_cols:
             print("Error: No sample columns found in TPM data after merging.")
             return None

        # Calculate sum of TPMs per gene for each sample
        # Group by gene, sum TPMs across all transcripts belonging to that gene
        gene_tpm_sum = merged_df.groupby('gene')[sample_cols].sum()

        # Add gene TPM sum back to the merged dataframe
        merged_df = pd.merge(merged_df, gene_tpm_sum, on='gene', how='left', suffixes=('', '_gene_sum'))

        # Calculate PSI for each transcript in each sample
        psi_df = merged_df[['transcript'] + sample_cols].copy() # Start with transcript id and sample tpm columns
        psi_df = psi_df.set_index('transcript') # Set transcript as index now

        for sample in sample_cols:
            gene_sum_col = sample + '_gene_sum'
            # PSI = transcript TPM / gene TPM sum
            # Avoid division by zero: if gene sum is 0, PSI is 0 or NaN (SUPPA uses 0)
            psi_df[sample] = merged_df.apply(
                lambda row: row[sample] / row[gene_sum_col] if row[gene_sum_col] > 0 else 0.0,
                axis=1
            ).values # Assign values based on original merged_df row alignment

        # Ensure the final PSI dataframe has transcripts as index and samples as columns
        # This should already be the case, but double-check formatting if needed.

        # Filter PSI df to only include transcripts originally in tpm_df index
        # psi_df = psi_df.loc[tpm_df.index.intersection(psi_df.index)] # Redundant if merge was inner

        return psi_df

    except Exception as e:
        print(f"Error during PSI calculation: {e}")
        import traceback
        traceback.print_exc()
        return None

def diffSplice_transcript_direct(psi_a_df, psi_b_df):
    """
    Calculates differential transcript usage (dPSI) per transcript between two conditions.

    Args:
        psi_a_df (pd.DataFrame): PSI values for condition A (transcripts x samples).
        psi_b_df (pd.DataFrame): PSI values for condition B (transcripts x samples).

    Returns:
        pd.DataFrame: Dataframe with transcript, dPSI (PSI_B - PSI_A), and p-value (Wilcoxon).
                      Returns None if calculation fails.
    """
    try:
        # Ensure indices are aligned (common transcripts)
        common_transcripts = psi_a_df.index.intersection(psi_b_df.index)
        if len(common_transcripts) == 0:
            print("Error: No common transcripts between PSI datasets A and B.")
            return None

        psi_a_filtered = psi_a_df.loc[common_transcripts]
        psi_b_filtered = psi_b_df.loc[common_transcripts]

        results = []

        # Check for sufficient samples for statistical test
        min_samples = 2 # Minimum samples per group for Wilcoxon often suggested
        if psi_a_filtered.shape[1] < min_samples or psi_b_filtered.shape[1] < min_samples:
            print(f"Warning: Fewer than {min_samples} samples in one or both conditions. P-values may be unreliable or not computed.")

        for transcript in common_transcripts:
            psi_a_values = psi_a_filtered.loc[transcript].values.flatten()
            psi_b_values = psi_b_filtered.loc[transcript].values.flatten()

            # Handle potential NaN values if necessary (e.g., remove them)
            psi_a_values = psi_a_values[~np.isnan(psi_a_values)]
            psi_b_values = psi_b_values[~np.isnan(psi_b_values)]

            # Calculate mean PSI for dPSI calculation
            mean_psi_a = np.mean(psi_a_values) if len(psi_a_values) > 0 else np.nan
            mean_psi_b = np.mean(psi_b_values) if len(psi_b_values) > 0 else np.nan

            # Calculate dPSI = mean(PSI_B) - mean(PSI_A)
            dpsi = mean_psi_b - mean_psi_a if not (np.isnan(mean_psi_a) or np.isnan(mean_psi_b)) else np.nan

            # Perform Wilcoxon rank-sum test if enough data points
            p_value = np.nan
            if len(psi_a_values) >= min_samples and len(psi_b_values) >= min_samples:
                 # Check for variance before test
                 if np.var(psi_a_values) > 1e-9 or np.var(psi_b_values) > 1e-9: # Avoid test if all values in a group are identical
                     try:
                         stat, p_value = mannwhitneyu(psi_a_values, psi_b_values, alternative='two-sided', nan_policy='omit')
                     except ValueError as ve:
                         # Handle cases like identical data where Wilcoxon fails
                         # print(f"Could not compute Wilcoxon for {transcript}: {ve}")
                         p_value = 1.0 # Or assign NaN, depending on how you want to treat this case
                 else:
                     # If no variance in one or both groups, can't do ranksums meaningfully
                     p_value = 1.0 # Or np.nan

            elif len(psi_a_values) > 0 and len(psi_b_values) > 0:
                # If fewer than min_samples but still some data, maybe return dPSI without p-value
                # Or indicate p-value is unreliable. Here we return NaN p-value.
                pass # p_value remains np.nan

            results.append({'transcript': transcript, 'dPSI': dpsi, 'p_value': p_value})

        if not results:
            print("Error: No results generated.")
            return None

        diff_results_df = pd.DataFrame(results)
        return diff_results_df.set_index('transcript')

    except Exception as e:
        print(f"Error during differential splicing calculation: {e}")
        import traceback
        traceback.print_exc()
        return None


