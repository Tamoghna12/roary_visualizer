import os
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
import argparse
import warnings
import numpy as np
from tqdm import tqdm

# Suppress FutureWarning
warnings.filterwarnings('ignore', category=FutureWarning)

def create_interactive_bar(df, output_path, folder_name):
    fig = px.bar(
        x=df.iloc[:, 0],
        y=df.iloc[:, 2],
        color=df.iloc[:, 0],
        title=f'Gene Category Distribution - {folder_name}'
    )
    fig.update_layout(
        xaxis_title='Gene Category',
        yaxis_title='Number of Genes',
        width=800,
        height=600
    )
    fig.write_html(output_path)

def create_static_bar(df, output_path, folder_name):
    sns.set_theme(style="whitegrid", palette="husl")
    plt.figure(figsize=(12, 6))
    plt.bar(df.iloc[:, 0], df.iloc[:, 1], color=sns.color_palette("husl", len(df)))
    plt.title(f'Gene Category Distribution - {folder_name}', fontsize=12, fontweight='bold')
    plt.xlabel('Gene Category', fontsize=10)
    plt.ylabel('Number of Genes', fontsize=10)
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_interactive_pie(roary, output_path, folder_name):
    total_strains = roary.shape[1]
    core = roary[(roary.sum(axis=1) >= total_strains * 0.99)].shape[0]
    softcore = roary[(roary.sum(axis=1) >= total_strains * 0.95) & (roary.sum(axis=1) < total_strains * 0.99)].shape[0]
    shell = roary[(roary.sum(axis=1) >= total_strains * 0.15) & (roary.sum(axis=1) < total_strains * 0.95)].shape[0]
    cloud = roary[roary.sum(axis=1) < total_strains * 0.15].shape[0]

    fig = go.Figure(data=[go.Pie(
        labels=['Core', 'Soft-core', 'Shell', 'Cloud'],
        values=[core, softcore, shell, cloud],
        textinfo='label+percent',
        insidetextorientation='radial'
    )])
    fig.update_layout(
        title=f'Pangenome Distribution - {folder_name}',
        width=800,
        height=800
    )
    fig.write_html(output_path)

def create_static_pie(roary, output_path, folder_name):
    total_strains = roary.shape[1]
    core = roary[(roary.sum(axis=1) >= total_strains * 0.99)].shape[0]
    softcore = roary[(roary.sum(axis=1) >= total_strains * 0.95) & (roary.sum(axis=1) < total_strains * 0.99)].shape[0]
    shell = roary[(roary.sum(axis=1) >= total_strains * 0.15) & (roary.sum(axis=1) < total_strains * 0.95)].shape[0]
    cloud = roary[roary.sum(axis=1) < total_strains * 0.15].shape[0]

    plt.figure(figsize=(10, 10))
    plt.pie(
        [core, softcore, shell, cloud],
        labels=['Core', 'Soft-core', 'Shell', 'Cloud'],
        autopct='%1.1f%%',
        startangle=90,
        colors=sns.color_palette("husl", 4)
    )
    plt.title(f'Pangenome Distribution - {folder_name}', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_static_matrix(roary, tree_file, output_path, folder_name, labels=False):
    try:
        t = Phylo.read(tree_file, 'newick')
    except Exception as e:
        print(f"Skipping {folder_name}: Error reading tree file - {str(e)}")
        return

    plt.figure(figsize=(17, 10))
    ax1 = plt.subplot2grid((1, 40), (0, 10), colspan=30)
    ax1.matshow(roary.T, cmap=plt.cm.Blues, vmin=0, vmax=1, aspect='auto', interpolation='none')
    ax1.axis('off')
    ax = plt.subplot2grid((1, 40), (0, 0), colspan=10, facecolor='white')
    plt.title(f'Tree\n({roary.shape[1]} strains)', fontsize=12, fontweight='bold')
    if labels:
        Phylo.draw(t, axes=ax, show_confidence=False, label_func=lambda x: str(x)[:10], do_show=False)
    else:
        Phylo.draw(t, axes=ax, show_confidence=False, label_func=lambda x: None, do_show=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def process_folder(folder_path, folder_name):
    summary_file = os.path.join(folder_path, 'summary_statistics.txt')
    tree_file = os.path.join(folder_path, 'accessory_binary_genes.fa.newick')
    spreadsheet_file = os.path.join(folder_path, 'gene_presence_absence.csv')

    if not os.path.exists(summary_file) or not os.path.exists(spreadsheet_file):
        print(f"Skipping {folder_name}: Missing required files")
        return

    try:
        df = pd.read_csv(summary_file, sep='\t', header=None)
        create_interactive_bar(df, os.path.join(folder_path, f'{folder_name}_barplot.html'), folder_name)
        create_static_bar(df, os.path.join(folder_path, f'{folder_name}_barplot.png'), folder_name)

        roary = pd.read_csv(spreadsheet_file, low_memory=False)
        roary.set_index('Gene', inplace=True)
        roary.drop(columns=roary.columns[:14], inplace=True)
        roary.replace('.{2,100}', 1, regex=True, inplace=True)
        roary.fillna(0, inplace=True)
        roary = roary.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

        create_interactive_pie(roary, os.path.join(folder_path, f'{folder_name}_pie.html'), folder_name)
        create_static_pie(roary, os.path.join(folder_path, f'{folder_name}_pie.png'), folder_name)

        if os.path.exists(tree_file):
            create_static_matrix(roary, tree_file, os.path.join(folder_path, f'{folder_name}_matrix.png'), folder_name, labels=True)
        else:
            print(f"Skipping matrix plots for {folder_name}: Missing tree file")

        print(f"Successfully processed {folder_name}")
    except Exception as e:
        print(f"Error processing {folder_name}: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description='Automate Roary plotting process')
    parser.add_argument('--input_dir', required=False, help='Path to directory containing Roary output folders')
    args = parser.parse_args()

    if not args.input_dir:
        input_directory = input("Enter the path to the Roary output directory: ").strip()
    else:
        input_directory = args.input_dir

    if not os.path.exists(input_directory) or not os.path.isdir(input_directory):
        print("Invalid directory path. Please check and try again.")
        return

    folders = [f for f in os.listdir(input_directory) if os.path.isdir(os.path.join(input_directory, f)) and '_' in f]
    for folder_name in tqdm(folders, desc="Processing folders"):
        folder_path = os.path.join(input_directory, folder_name)
        process_folder(folder_path, folder_name)

    print("\nProcessing complete!")

if __name__ == "__main__":
    main()
