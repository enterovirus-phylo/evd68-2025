# A cleaned up R file that plots a few different plots. Only the last is used in the manuscript.
# These generally plot the distribution of subgenotype and number of sequences - either by year
# or by season 

# Produces: 
# 1.  A plot of the subgenogroup composition by season (winter/summer) from 2013-2025, by number of
#   sequences. This is *not* used in the manuscript.
# 2.  A stacked plot showing number of sequences per season and the subgenotype composition by season
#   This is also *not* used in the manuscript.
# 3.  A stacked bar plot showing number of sequences per calendar year and the subgenotype composition 
#   by year. This *IS* used in the manuscript.


import pandas as pd
import matplotlib.pyplot as plt

## Run nextclade on the 300 sequence build first, to get proper subgenotype assignments
#
# To do this Nextclade Run:
# From within wsl/enteroviruses
# nextclade run --input-dataset nextclade-datasets/ev-d68 \
#  --output-all nextclade-results/ev-d68-300bp \
#  ev_githubs/enterovirus_d68-oct2025/vp1/results/aligned_300.fasta

####
# To run THIS SCRIPT:
## run from enteroviruses/ev_githubs/enterovirus_d68-oct2025 directory
# like this:
# python scripts/plotSequencesByGenogroupSeason.py 
# --meta-in vp1/results/metadata_subgeno_seqlen_dupcol_extramet.tsv --meta300-in 
# vp1/results/metadata_subgeno_seqlen_dupcol_extramet_300-23Nov2025.tsv --nextclade-in  ../../nextclade-results/ev-d68-300bp/nextclade.tsv 
# --fig-out vp1/results/newfig.png

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='parse metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--meta-in', help="input final meta file - fall 2025 build")
    parser.add_argument('--meta300-in', help="run with >=300bp sequences")
    parser.add_argument('--nextclade-in', help="nextclade results tsv file for 300bp sequences")
    parser.add_argument('--fig-out', help="final figure file")
    args = parser.parse_args()
    
    #for ipython debugging
    #meta_in = "vp1/results/metadata_subgeno_seqlen_dupcol_extramet.tsv"
    #meta300_in = "vp1/results/metadata_subgeno_seqlen_dupcol_extramet_300-23Nov2025.tsv"
    #nextclade_in = "../../nextclade-results/ev-d68-300bp/nextclade.tsv"

    
    # read in metadata file
    meta = pd.read_csv(args.meta_in, sep='\t', index_col=False)
    # We want to filter the metadata to only European region
    meta2 = meta[meta['region'] == 'Europe']
     #Filter to only metadata with date after 2013-08-01
    meta2 = meta2[meta2['date'] >= '2013-08-01']
    # Note that some dates may be missing month and be in format 2022-XX-XX - we should exclude these first
    meta2 = meta2[~meta2['date'].str.contains('XX-XX')]

    # At least for Oct 2025 data - seq & metadata length of file files matches, so no need
    # to use one to check that only final data is in there - both are 3771 which matches
    # what's in the build
    
    # Metadata file has columns age (exact age), region, country, subgenogroup, date
    # and age_range1, age_range2, age_range3 (decreasing bin granularity, with age_range3
    # having the most data as it's least granular)

    # read in 300bp metadata file
    meta300 = pd.read_csv(args.meta300_in, sep='\t', index_col=False)
    # using seq-len column, extract all sequences that are < 700bp 
    meta300_short = meta300[meta300['seq-len'] < 700]
    # take only those from Europe (loose about 600)
    meta300_short = meta300_short[meta300_short['region'] == 'Europe']
    #Filter to only metadata with date after 2013-08-01
    meta300_short = meta300_short[meta300_short['date'] >= '2013-08-01']
    meta300_short = meta300_short[~meta300_short['date'].str.contains('XX-XX')]

    #count how many are genogroup B1 and how many are B3 in meta300_short based on the 'subgenogroup' column - to check that it looks reasonable
    print(meta300_short['subgenogroup'].value_counts())

    # We can't trust the genogroup assignments here, since B1 was incorrectly called on the tree
    # So let's go get it from the nextclade results and replace the subgenogroup column
    nextclade_results = pd.read_csv(args.nextclade_in, sep='\t')
    # find the rows in nextclade_results as 'seqName' that are also in meta300_short as 'strain'
    matching_rows = nextclade_results[nextclade_results['seqName'].isin(meta300_short['strain'])]
    # Now we can create a mapping from strain to subgenogroup using the 'clade' column in nextclade results
    strain_to_subgenogroup = dict(zip(matching_rows['seqName'], matching_rows['clade']))
    # count how many of these are B1 and how many are B3 to check that it looks reasonable
    print(matching_rows['clade'].value_counts())
    # Now we can replace the 'subgenogroup' column in meta300_short with the values from this mapping
    meta300_short.loc[:, 'subgenogroup'] = meta300_short['strain'].map(strain_to_subgenogroup)

    # Combine meta300_short and meta2 - to incorporate all VP1 sequences that are >=300bp
    # We can do this by concatenating the two dataframes and then dropping duplicates based on the 'accession' column
    meta_combined = pd.concat([meta2, meta300_short], ignore_index=True)
    meta_combined2 = meta_combined.drop_duplicates(subset='strain')

    #rename any subgenogroup 'A2/D' to be A2  
    meta_combined2.loc[meta_combined2['subgenogroup'] == 'A2/D', 'subgenogroup'] = 'A2'

    # We need seasons from February to July and August to January each year
    # So for example February to July 2014 and August 2014 to January 2015.
    # There may also be dates with missing day, in format 2022-05-XX - we can keep
    # these since we are sorting only by month, but need to make sure to handle them correctly when assigning seasons.
    # We can create a new column 'season' to indicate which season each sample belongs to
    # Handle dates with missing day in format 2022-04-XX  
    def assign_season(date):
        month = int(date.split('-')[1])
        year = int(date.split('-')[0])
        if month >= 2 and month <= 8:
            return f"{year}"
        elif month >= 9 and month <= 12:
            return f"{year}-{year+1}"
        else: # month is January
            return f"{year-1}-{year}"

    # Avoid SettingsWithCopyWarning by using .loc to assign the new column
    meta_combined2.loc[:, 'season'] = meta_combined2['date'].apply(assign_season)

    # Now we can group by season and count the number of samples in each season
    season_counts = meta_combined2['season'].value_counts().sort_index()

    # Now group by both season and subgenogroup to get counts of each subgenogroup in each season
    season_subgenogroup_counts = meta_combined2.groupby(['season', 'subgenogroup']).size().unstack(fill_value=0)

    # Rename seasons into 'Summer 2014' and 'Winter 2014' etc. for better plotting
    def rename_season(season): 
        if '-' in season: 
            year1, year2 = season.split('-')
            return f"Winter {year1}"
        else: 
            return f"Summer {season}"
    #apply this to season_subgenogroup_counts index
    season_subgenogroup_counts.index = season_subgenogroup_counts.index.map(rename_season)

    # Add a row for 'Winter 2020' with all counts 0 to ensure it appears in the plot even if there are no samples from that season
    # Ensure it gets added after Summer 2020 and before Summer 2021
    season_subgenogroup_counts.loc['Winter 2020'] = 0
    season_subgenogroup_counts = season_subgenogroup_counts.sort_index(key=lambda x: x.map(lambda s: (int(s.split()[1]), 0 if s.startswith('Summer') else 1)))

    # Now create a plot of this data - use side-by-side bar plots with season on the x-axis, counts on the y-axis, and different bars for each subgenogroup
    # Plot to a file specified by --fig-out
    # Ensure X-axis labels aren't getting cut off - angle them and put more room at the bottom of the plot
    season_subgenogroup_counts.plot(kind='bar', stacked=False, figsize=(10,6))
    plt.xlabel('Season')
    plt.ylabel('Number of Samples')
    plt.title('Number of Samples by Season and Subgenogroup in Europe')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig("figs/season_subgenogroup_counts.png")
    #plt.savefig(args.fig_out)
    plt.close()
    ### NOTE THIS PLOT NOT USED IN THE MANUSCRIPT

    ###############################################
    ###############################################
    # Two-panel figure showing total sequences per season and subgenotype composition
    import matplotlib.gridspec as gridspec
    import numpy as np

    # Create the complete list of all seasons from Winter 2013 to Winter 2025
    all_seasons_list = []
    for yr in range(2013, 2025):
        all_seasons_list.append(f"Winter {yr}")
        all_seasons_list.append(f"Summer {yr + 1}")
    all_seasons_list.append("Winter 2025")

    # Use the existing season_subgenogroup_counts dataframe
    # Add missing seasons with all zeros
    season_subgeno_counts = season_subgenogroup_counts.copy()
    for season in all_seasons_list:
        if season not in season_subgeno_counts.index:
            season_subgeno_counts.loc[season] = 0

    # Reindex to ensure correct chronological order
    season_subgeno_counts = season_subgeno_counts.reindex(all_seasons_list)

    # Ensure all subgenotypes are present as columns
    all_subgenotypes = ['A1', 'A2', 'B1', 'B2', 'B3']
    for subgeno in all_subgenotypes:
        if subgeno not in season_subgeno_counts.columns:
            season_subgeno_counts[subgeno] = 0

    # Select and reorder columns to match desired subgenotype order
    season_subgeno_counts = season_subgeno_counts[all_subgenotypes]

    # Calculate total counts per season
    season_counts_total = season_subgeno_counts.sum(axis=1)

    # Calculate proportions (handle seasons with zero counts)
    season_subgeno_props = season_subgeno_counts.div(season_counts_total, axis=0).fillna(0)

    # Define consistent colors for subgenotypes
    colors = {
        'A1': '#1f77b4',
        'A2': '#ff7f0e',
        'B1': '#2ca02c',
        'B2': '#d62728',
        'B3': '#9467bd'
    }
    color_list = [colors[subgeno] for subgeno in all_subgenotypes]

    # Create figure with two panels - sized for publication
    fig = plt.figure(figsize=(11, 5.0))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3], hspace=0.45)

    # Top panel: total counts per season
    ax_top = fig.add_subplot(gs[0])
    x_pos = np.arange(len(season_counts_total))
    ax_top.bar(x_pos, season_counts_total.values, color="#566357", width=0.75, edgecolor='black', linewidth=0.6)
    ax_top.set_ylabel('Number of\nSequences', fontsize=11)
    ax_top.set_title('Total Sequences per Season', fontsize=12, fontweight='bold')
    ax_top.set_xticks(x_pos)
    ax_top.set_xticklabels([])  # Empty labels; shared x-axis with bottom panel
    ax_top.tick_params(axis='y', labelsize=10)
    ax_top.grid(axis='y', alpha=0.3, linestyle='--')

    # Bottom panel: 100%-stacked bar chart of subgenotype proportions
    ax_bottom = fig.add_subplot(gs[1])
    bottom_vals = np.zeros(len(season_subgeno_props))

    for i, subgeno in enumerate(all_subgenotypes):
        ax_bottom.bar(x_pos, season_subgeno_props[subgeno].values,
                      bottom=bottom_vals, label=subgeno, color=color_list[i],
                      width=0.75, edgecolor='black', linewidth=0.6)
        bottom_vals += season_subgeno_props[subgeno].values

    ax_bottom.set_xlabel('Season', fontsize=11)
    ax_bottom.set_ylabel('Proportion', fontsize=11)
    ax_bottom.set_title('Subgenotype Composition by Season', fontsize=12, fontweight='bold')
    ax_bottom.set_xticks(x_pos)
    ax_bottom.set_xticklabels(all_seasons_list, rotation=45, ha='right', fontsize=10)
    ax_bottom.tick_params(axis='y', labelsize=10)
    ax_bottom.set_ylim(0, 1)
    ax_bottom.legend(title='Subgenotype', loc='upper left', bbox_to_anchor=(1.01, 1),
                     fontsize=10, title_fontsize=11, frameon=True)
    ax_bottom.grid(axis='y', alpha=0.3, linestyle='--')

    plt.tight_layout(pad=1.0)
    plt.subplots_adjust(right=0.85)  # Remove extra space on right side
    #plt.savefig("vp1/results/season_composition_two_panel.png", dpi=300, bbox_inches='tight')
    plt.savefig("figs/season_composition_two_panel.png", dpi=300, bbox_inches='tight')
    plt.close()




    ###############################################
    ###############################################
    # Two-panel figure showing total sequences per calendar year and subgenotype composition

    # Extract year from date column (handles dates with missing days like 2022-04-XX)
    meta_combined2['year'] = meta_combined2['date'].str[:4].astype(int)

    # Get all years in the dataset (from 2013 onwards)
    min_year = meta_combined2['year'].min()
    max_year = meta_combined2['year'].max()
    all_years_list = list(range(int(min_year), int(max_year) + 1))

    # Check subgenogroup values - print unique values and count missing
    print("\nSubgenogroup value counts:")
    print(meta_combined2['subgenogroup'].value_counts(dropna=False))
    print(f"\nRows with missing subgenogroup: {meta_combined2['subgenogroup'].isna().sum()}")

    # Filter to only rows with valid subgenotypes (A1, A2, B1, B2, B3)
    valid_subgenotypes = ['A1', 'A2', 'B1', 'B2', 'B3']
    meta_combined2_filtered = meta_combined2[meta_combined2['subgenogroup'].isin(valid_subgenotypes)]
    print(f"Rows with valid subgenotypes: {len(meta_combined2_filtered)} out of {len(meta_combined2)}")
    print(meta_combined2[~meta_combined2['subgenogroup'].isin(['A1', 'A2', 'B1', 'B2', 'B3'])]['subgenogroup'].value_counts())

    # Group by year to get total counts per year (using filtered data)
    year_counts_total = meta_combined2_filtered.groupby('year').size()

    # Group by year and subgenogroup
    year_subgeno_counts = meta_combined2_filtered.groupby(['year', 'subgenogroup']).size().unstack(fill_value=0)

    # Ensure all years are present (fill missing years with 0)
    year_subgeno_counts = year_subgeno_counts.reindex(all_years_list, fill_value=0)
    year_counts_total = year_counts_total.reindex(all_years_list, fill_value=0)

    # Ensure all subgenotypes are present as columns
    for subgeno in all_subgenotypes:
        if subgeno not in year_subgeno_counts.columns:
            year_subgeno_counts[subgeno] = 0

    # Select and reorder columns to match desired subgenotype order
    year_subgeno_counts = year_subgeno_counts[all_subgenotypes]
    #rename A2 to A2/D if it exists in the data - to match the colors defined below
    if 'A2' in year_subgeno_counts.columns:
        year_subgeno_counts.rename(columns={'A2': 'A2/D'}, inplace=True)
        all_subgenotypes = ['A1', 'A2/D', 'B1', 'B2', 'B3']

    # Calculate proportions (handle years with zero counts)
    year_subgeno_props = year_subgeno_counts.div(year_counts_total, axis=0).fillna(0)

    # Define new consistent colors for subgenotypes
    colors = {
        'A1': '#5fb3f3',
        'A2/D': '#7fdbbd',
        'B1': '#e8ed5d',
        'B2': '#deae3d',
        'B3': '#ff7e0c'
    }
    color_list = [colors[subgeno] for subgeno in all_subgenotypes]

    # Create figure with two panels - same sizing as seasonal plot
    fig = plt.figure(figsize=(11, 5.0))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3], hspace=0.45)

    # Top panel: total counts per year
    ax_top = fig.add_subplot(gs[0])
    x_pos = np.arange(len(year_counts_total))
    ax_top.bar(x_pos, year_counts_total.values, color="#566357", width=0.75, edgecolor='black', linewidth=0.6)
    ax_top.set_ylabel('Number of\nSequences', fontsize=11)
    ax_top.set_title('Total Sequences per Year', fontsize=12, fontweight='bold')
    ax_top.set_xticks(x_pos)
    ax_top.set_xticklabels([])  # Empty labels; shared x-axis with bottom panel
    ax_top.tick_params(axis='y', labelsize=10)
    ax_top.grid(axis='y', alpha=0.3, linestyle='--')

    # Bottom panel: 100%-stacked bar chart of subgenotype proportions
    ax_bottom = fig.add_subplot(gs[1])
    bottom_vals = np.zeros(len(year_subgeno_props))

    for i, subgeno in enumerate(all_subgenotypes):
        ax_bottom.bar(x_pos, year_subgeno_props[subgeno].values,
                      bottom=bottom_vals, label=subgeno, color=color_list[i],
                      width=0.75, edgecolor='black', linewidth=0.6)
        bottom_vals += year_subgeno_props[subgeno].values

    ax_bottom.set_xlabel('Year', fontsize=11)
    ax_bottom.set_ylabel('Proportion', fontsize=11)
    ax_bottom.set_title('Clade Composition by Year', fontsize=12, fontweight='bold')
    ax_bottom.set_xticks(x_pos)
    ax_bottom.set_xticklabels(all_years_list, rotation=45, ha='right', fontsize=10)
    ax_bottom.tick_params(axis='y', labelsize=10)
    ax_bottom.set_ylim(0, 1)
    ax_bottom.legend(title='Clade', loc='upper left', bbox_to_anchor=(1.01, 1),
                     fontsize=10, title_fontsize=11, frameon=True)
    ax_bottom.grid(axis='y', alpha=0.3, linestyle='--')

    plt.tight_layout(pad=1.0)
    plt.subplots_adjust(right=0.85)  # Remove extra space on right side
    #plt.savefig("vp1/results/year_composition_two_panel.png", dpi=300, bbox_inches='tight')
    plt.savefig("figs/year_composition_two_panel.png", dpi=300, bbox_inches='tight')
    plt.close()

 