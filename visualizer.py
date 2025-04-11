import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from Bio import Phylo
import os
import tempfile
import shutil
import io
import base64
from typing import Dict, Any, Optional, Tuple
import re

# Configure the page
st.set_page_config(
    page_title="Roary Visualizer",
    layout="wide",
    page_icon="🧬"
)

# Initialize session state for persistent variables
if 'temp_dir' not in st.session_state:
    st.session_state.temp_dir = tempfile.mkdtemp()

# Function to validate file formats
def validate_file(file, expected_extension):
    """Validate if a file has the expected extension
    
    Args:
        file: The uploaded file object
        expected_extension: String with expected file extension
        
    Returns:
        bool: True if file is valid, False otherwise
    """
    filename = file.name.lower()
    if not filename.endswith(expected_extension.lower()):
        return False
    return True

# File processing functions
@st.cache_data
def load_summary_statistics(file_path: str) -> Optional[pd.DataFrame]:
    """Load and parse summary statistics file
    
    Args:
        file_path: Path to summary_statistics.txt
        
    Returns:
        DataFrame with parsed statistics or None if error
    """
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
        if df.shape[1] < 3:  # Basic validation
            st.error("Summary statistics file appears to be malformed.")
            return None
        return df
    except Exception as e:
        st.error(f"Error loading summary statistics: {str(e)}")
        return None

@st.cache_data
def load_gene_presence_absence(file_path: str) -> Optional[pd.DataFrame]:
    """Load and process gene presence/absence matrix
    
    Args:
        file_path: Path to gene_presence_absence.csv
        
    Returns:
        Processed DataFrame or None if error
    """
    try:
        roary = pd.read_csv(file_path, low_memory=False)
        
        # Validate that required columns exist
        if 'Gene' not in roary.columns:
            st.error("Gene column not found in gene_presence_absence.csv")
            return None
            
        # Set gene as index and drop annotation columns
        roary.set_index('Gene', inplace=True)
        
        # Find first genome column (assuming annotation columns come before genomes)
        genome_col_start = 0
        for i, col in enumerate(roary.columns):
            if re.match(r'^\w+\.\w+$', col):  # Simple pattern for genome identifiers
                genome_col_start = i
                break
        
        if genome_col_start > 0:
            roary.drop(list(roary.columns[:genome_col_start]), axis=1, inplace=True)
        else:
            # If no genome columns are identified, use old method but warn the user
            st.warning("Could not clearly identify genome columns, using first 14 columns as annotation.")
            roary.drop(list(roary.columns[:14]), axis=1, inplace=True)
        
        # Convert presence/absence data
        # Replace any non-zero strings with 1 (presence)
        roary = roary.applymap(lambda x: 1 if isinstance(x, str) and len(x) > 0 else x)
        roary.fillna(0, inplace=True)
        
        # Convert to numeric, more efficient than looping
        for col in roary.columns:
            roary[col] = pd.to_numeric(roary[col], errors='coerce').fillna(0).astype(int)
        
        return roary
    except Exception as e:
        st.error(f"Error loading gene presence/absence matrix: {str(e)}")
        return None

@st.cache_data
def load_newick_tree(file_path: str) -> Optional[Any]:
    """Load phylogenetic tree from newick file
    
    Args:
        file_path: Path to newick tree file
        
    Returns:
        Loaded tree object or None if error
    """
    if not os.path.exists(file_path):
        return None
        
    try:
        tree = Phylo.read(file_path, 'newick')
        return tree
    except Exception as e:
        st.error(f"Error loading phylogenetic tree: {str(e)}")
        return None

def process_roary_files(file_paths: Dict[str, str]) -> Dict[str, Any]:
    """Process all Roary output files
    
    Args:
        file_paths: Dictionary of file names to file paths
        
    Returns:
        Dictionary containing processed data objects
    """
    results = {}
    
    # Load summary statistics
    if 'summary_statistics.txt' in file_paths:
        results['df'] = load_summary_statistics(file_paths['summary_statistics.txt'])
    else:
        st.warning("Summary statistics file (summary_statistics.txt) not found.")
    
    # Load gene presence/absence
    if 'gene_presence_absence.csv' in file_paths:
        results['roary'] = load_gene_presence_absence(file_paths['gene_presence_absence.csv'])
    else:
        st.warning("Gene presence/absence file (gene_presence_absence.csv) not found.")
    
    # Load tree if available
    if 'accessory_binary_genes.fa.newick' in file_paths:
        results['tree'] = load_newick_tree(file_paths['accessory_binary_genes.fa.newick'])
        results['tree_file'] = file_paths['accessory_binary_genes.fa.newick']
    
    return results

# Add download buttons for various file formats
def download_buttons(fig, base_filename: str):
    """Create download buttons for PNG and SVG formats
    
    Args:
        fig: Plotly or Matplotlib figure object
        base_filename: Base filename without extension
    """
    col1, col2 = st.columns(2)
    
    # PNG Download
    with col1:
        png_buffer = io.BytesIO()
        if 'matplotlib.figure.Figure' in str(type(fig)):
            fig.savefig(png_buffer, format="png", dpi=300, bbox_inches='tight')
        else:  # Plotly figure
            fig.write_image(png_buffer, format="png")
        png_buffer.seek(0)
        st.download_button(
            label="Download as PNG",
            data=png_buffer,
            file_name=f"{base_filename}.png",
            mime="image/png"
        )
    
    # SVG Download
    with col2:
        svg_buffer = io.BytesIO()
        if 'matplotlib.figure.Figure' in str(type(fig)):
            fig.savefig(svg_buffer, format="svg", bbox_inches='tight')
        else:  # Plotly figure
            fig.write_image(svg_buffer, format="svg")
        svg_buffer.seek(0)
        st.download_button(
            label="Download as SVG",
            data=svg_buffer,
            file_name=f"{base_filename}.svg",
            mime="image/svg+xml"
        )

# Visualization functions
def plot_bar(df: pd.DataFrame, log_scale: bool = False) -> None:
    """Create a bar plot of gene categories
    
    Args:
        df: DataFrame with gene categories
        log_scale: Whether to use logarithmic scale for y-axis
    """
    try:
        fig = px.bar(
            x=df.iloc[:, 0],
            y=df.iloc[:, 2],
            color=df.iloc[:, 0],
            labels={'x': 'Gene Category', 'y': 'Number of Genes'},
            title='Gene Category Distribution'
        )
        
        if log_scale:
            fig.update_layout(yaxis_type="log")
            
        fig.update_layout(
            xaxis_title="Gene Category",
            yaxis_title="Number of Genes",
            legend_title="Gene Category",
            height=500
        )
        
        # Display the figure
        st.plotly_chart(fig, use_container_width=True)
        
        # Add download buttons
        download_buttons(fig, "gene_category_distribution")
        
    except Exception as e:
        st.error(f"Error creating bar plot: {str(e)}")

def plot_pie(roary: pd.DataFrame, core_threshold: float, softcore_threshold: float, shell_threshold: float) -> None:
    """Create a pie chart of pangenome distribution
    
    Args:
        roary: Processed gene presence/absence matrix
        core_threshold: Threshold for core genes (proportion of genomes)
        softcore_threshold: Threshold for softcore genes
        shell_threshold: Threshold for shell genes
    """
    try:
        total_strains = roary.shape[1]
        
        core = roary[(roary.sum(axis=1) >= total_strains * core_threshold)].shape[0]
        softcore = roary[(roary.sum(axis=1) >= total_strains * softcore_threshold) & 
                         (roary.sum(axis=1) < total_strains * core_threshold)].shape[0]
        shell = roary[(roary.sum(axis=1) >= total_strains * shell_threshold) & 
                      (roary.sum(axis=1) < total_strains * softcore_threshold)].shape[0]
        cloud = roary[roary.sum(axis=1) < total_strains * shell_threshold].shape[0]
        
        fig = go.Figure(data=[go.Pie(
            labels=['Core', 'Soft-core', 'Shell', 'Cloud'],
            values=[core, softcore, shell, cloud],
            textinfo='label+percent',
            marker=dict(colors=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
        )])
        
        fig.update_layout(
            title={'text': 'Pangenome Distribution', 'y':0.92},
            height=500
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Add download buttons
        download_buttons(fig, "pangenome_distribution")
        
        # Display counts in a table
        st.markdown("### Gene Count by Category")
        count_data = {
            "Category": ["Core", "Soft-core", "Shell", "Cloud", "Total"],
            "Count": [core, softcore, shell, cloud, core+softcore+shell+cloud],
            "Percentage": [
                f"{core/(core+softcore+shell+cloud)*100:.1f}%",
                f"{softcore/(core+softcore+shell+cloud)*100:.1f}%",
                f"{shell/(core+softcore+shell+cloud)*100:.1f}%",
                f"{cloud/(core+softcore+shell+cloud)*100:.1f}%",
                "100%"
            ]
        }
        st.table(pd.DataFrame(count_data))
        
        # Add option to download the table data
        csv_buffer = io.StringIO()
        pd.DataFrame(count_data).to_csv(csv_buffer, index=False)
        st.download_button(
            label="Download Table Data (CSV)",
            data=csv_buffer.getvalue(),
            file_name="pangenome_category_counts.csv",
            mime="text/csv"
        )
        
    except Exception as e:
        st.error(f"Error creating pie chart: {str(e)}")

def plot_matrix(roary: pd.DataFrame, tree: Any, cluster_rows: bool = True, cluster_cols: bool = False) -> None:
    """Create a presence/absence matrix heatmap
    
    Args:
        roary: Processed gene presence/absence matrix
        tree: Phylogenetic tree object
        cluster_rows: Whether to cluster rows
        cluster_cols: Whether to cluster columns
    """
    try:
        plt.rcParams['svg.fonttype'] = 'none'  # Ensure text is editable in SVG
        
        if cluster_rows and cluster_cols:
            # For clustered heatmap, use clustermap instead of heatmap
            g = sns.clustermap(roary.T, cmap='Blues', figsize=(15, 10))
            fig = g.fig
            fig.suptitle("Presence/Absence Matrix (Clustered)", fontsize=16)
            plt.tight_layout()
        else:
            fig, ax = plt.subplots(figsize=(15, 10))
            
            if tree:
                # If we have a tree, use it to sort the columns
                leaf_names = [tip.name for tip in tree.get_terminals()]
                # Filter to only include names present in the dataframe
                valid_names = [name for name in leaf_names if name in roary.columns]
                
                if valid_names:
                    roary_sorted = roary[valid_names].T
                else:
                    roary_sorted = roary.T
                    st.warning("Tree leaf names do not match dataset genome names. Using original order.")
            else:
                roary_sorted = roary.T
                
            sns.heatmap(roary_sorted, cmap='Blues', cbar=True, ax=ax)
            ax.set_title("Presence/Absence Matrix")
            ax.set_xlabel("Genes")
            ax.set_ylabel("Genomes")
            plt.tight_layout()
        
        st.pyplot(fig)
        
        # Add download buttons
        download_buttons(fig, "presence_absence_matrix")
        
    except Exception as e:
        st.error(f"Error creating matrix plot: {str(e)}")

def generate_frequency_plot(roary: pd.DataFrame) -> None:
    """Generate a gene frequency plot
    
    Args:
        roary: Processed gene presence/absence matrix
    """
    try:
        # Calculate gene frequency (how many strains have each gene)
        gene_freq = roary.sum(axis=1)
        freq_count = gene_freq.value_counts().sort_index()
        
        fig = px.bar(
            x=freq_count.index,
            y=freq_count.values,
            labels={'x': 'Number of Genomes', 'y': 'Number of Genes'},
            title='Gene Frequency Plot'
        )
        
        fig.update_layout(
            xaxis_title="Number of Genomes",
            yaxis_title="Number of Genes",
            height=500
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Add download buttons
        download_buttons(fig, "gene_frequency_plot")
        
        # Add option to display and download frequency data
        if st.checkbox("Show frequency data table"):
            # Create a dataframe of the frequency data
            freq_df = pd.DataFrame({
                'Number_of_Genomes': freq_count.index,
                'Number_of_Genes': freq_count.values
            })
            st.dataframe(freq_df)
            
            # Download button for frequency data
            csv_buffer = io.StringIO()
            freq_df.to_csv(csv_buffer, index=False)
            st.download_button(
                label="Download Frequency Data (CSV)",
                data=csv_buffer.getvalue(),
                file_name="gene_frequency_data.csv",
                mime="text/csv"
            )
        
    except Exception as e:
        st.error(f"Error creating frequency plot: {str(e)}")

def generate_rarefaction_curve(roary: pd.DataFrame, num_permutations: int = 100) -> None:
    """Generate a rarefaction curve (gene accumulation curve)
    
    Args:
        roary: Processed gene presence/absence matrix
        num_permutations: Number of random permutations to generate
    """
    try:
        total_genomes = roary.shape[1]
        total_genes = roary.shape[0]
        
        # Limit permutations for large datasets
        if total_genomes > 50:
            num_permutations = min(num_permutations, 20)
            st.info(f"Large dataset detected. Limiting to {num_permutations} permutations for performance.")
        
        # Calculate rarefaction curve with multiple permutations
        genome_counts = list(range(1, total_genomes + 1))
        all_curves = []
        
        with st.spinner(f"Generating rarefaction curve with {num_permutations} permutations..."):
            for _ in range(num_permutations):
                # Shuffle genome order
                shuffled_cols = np.random.permutation(roary.columns)
                gene_counts = []
                
                # Calculate cumulative gene count
                gene_set = set()
                for i in range(total_genomes):
                    genome = shuffled_cols[i]
                    # Add genes present in this genome
                    present_genes = roary.index[roary[genome] > 0].tolist()
                    gene_set.update(present_genes)
                    gene_counts.append(len(gene_set))
                
                all_curves.append(gene_counts)
        
        # Calculate mean and std
        curve_array = np.array(all_curves)
        mean_curve = np.mean(curve_array, axis=0)
        std_curve = np.std(curve_array, axis=0)
        
        # Create plot
        fig = go.Figure()
        
        # Add mean line
        fig.add_trace(go.Scatter(
            x=genome_counts,
            y=mean_curve,
            mode='lines',
            name='Average',
            line=dict(color='blue', width=2)
        ))
        
        # Add confidence interval
        fig.add_trace(go.Scatter(
            x=genome_counts + genome_counts[::-1],
            y=np.concatenate([mean_curve + std_curve, (mean_curve - std_curve)[::-1]]),
            fill='toself',
            fillcolor='rgba(0,0,255,0.2)',
            line=dict(color='rgba(255,255,255,0)'),
            name='±1 Std Dev'
        ))
        
        fig.update_layout(
            title='Gene Accumulation Curve (Rarefaction)',
            xaxis_title='Number of Genomes',
            yaxis_title='Number of Genes',
            height=500
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Add download buttons
        download_buttons(fig, "rarefaction_curve")
        
        # Add rarefaction data download
        rarefaction_df = pd.DataFrame({
            'Genomes': genome_counts,
            'Mean_Genes': mean_curve,
            'Std_Dev': std_curve
        })
        
        csv_buffer = io.StringIO()
        rarefaction_df.to_csv(csv_buffer, index=False)
        st.download_button(
            label="Download Rarefaction Data (CSV)",
            data=csv_buffer.getvalue(),
            file_name="rarefaction_data.csv",
            mime="text/csv"
        )
        
    except Exception as e:
        st.error(f"Error creating rarefaction curve: {str(e)}")

def export_data(roary: pd.DataFrame) -> None:
    """Export processed data as CSV
    
    Args:
        roary: Processed gene presence/absence matrix
    """
    try:
        st.markdown("#### Download Options")
        
        # Create binary presence/absence matrix CSV
        csv_buffer = io.StringIO()
        roary.to_csv(csv_buffer)
        
        st.download_button(
            label="Download Processed Presence/Absence Matrix (CSV)",
            data=csv_buffer.getvalue(),
            file_name="processed_presence_absence.csv",
            mime="text/csv"
        )
        
        # Create gene frequency data
        gene_freq = pd.DataFrame({
            'Gene': roary.index,
            'Frequency': roary.sum(axis=1),
            'Present_in_genomes': roary.sum(axis=1) / roary.shape[1] * 100
        })
        
        freq_buffer = io.StringIO()
        gene_freq.to_csv(freq_buffer, index=False)
        
        st.download_button(
            label="Download Gene Frequency Data (CSV)",
            data=freq_buffer.getvalue(),
            file_name="gene_frequency.csv",
            mime="text/csv"
        )
        
        # Add binary matrix in TSV format
        tsv_buffer = io.StringIO()
        roary.to_csv(tsv_buffer, sep='\t')
        
        st.download_button(
            label="Download Presence/Absence Matrix (TSV)",
            data=tsv_buffer.getvalue(),
            file_name="presence_absence_matrix.tsv",
            mime="text/tab-separated-values"
        )
        
    except Exception as e:
        st.error(f"Error exporting data: {str(e)}")

# Main application layout
def main():
    """Main application function"""
    
    st.title("🧬 Roary Output Visualizer")
    st.markdown("""
    This application visualizes the output files from [Roary](https://sanger-pathogens.github.io/Roary/), 
    a tool for pan-genome analysis. Upload your Roary output files to generate interactive visualizations.
    """)
    
    # File upload section with better guidance
    st.markdown("### Step 1: Upload Your Roary Output Files")
    st.markdown("""
    Please upload the following files from your Roary output:
    - **summary_statistics.txt** (required): Contains gene category counts
    - **gene_presence_absence.csv** (required): Contains the gene presence/absence matrix
    - **accessory_binary_genes.fa.newick** (optional): Phylogenetic tree for strain ordering
    """)
    
    uploaded_files = st.file_uploader(
        "Upload Roary output files",
        accept_multiple_files=True,
        type=["txt", "csv", "newick"]
    )
    
    # Process uploaded files
    file_paths = {}
    if uploaded_files:
        for uploaded_file in uploaded_files:
            # Validate expected file extensions
            if uploaded_file.name == 'summary_statistics.txt' and not validate_file(uploaded_file, '.txt'):
                st.error(f"File {uploaded_file.name} doesn't appear to be a valid text file.")
                continue
                
            if uploaded_file.name == 'gene_presence_absence.csv' and not validate_file(uploaded_file, '.csv'):
                st.error(f"File {uploaded_file.name} doesn't appear to be a valid CSV file.")
                continue
                
            if 'newick' in uploaded_file.name and not validate_file(uploaded_file, '.newick'):
                st.error(f"File {uploaded_file.name} doesn't appear to be a valid Newick file.")
                continue
            
            path = os.path.join(st.session_state.temp_dir, uploaded_file.name)
            file_paths[uploaded_file.name] = path
            
            with open(path, 'wb') as f:
                f.write(uploaded_file.read())
    
    # Check if required files are present
    required_files = ['summary_statistics.txt', 'gene_presence_absence.csv']
    missing_files = [f for f in required_files if f not in file_paths]
    
    if missing_files:
        st.warning(f"Missing required files: {', '.join(missing_files)}")
        
    # Process files if required files are present
    if not missing_files and file_paths:
        results = process_roary_files(file_paths)
        
        if 'df' in results and results['df'] is not None and 'roary' in results and results['roary'] is not None:
            st.success("✅ Files loaded successfully!")
            
            # Basic dataset information
            st.markdown("### Dataset Information")
            st.markdown(f"- **Number of genomes**: {results['roary'].shape[1]}")
            st.markdown(f"- **Total genes**: {results['roary'].shape[0]}")
            
            # Visualization settings
            st.markdown("### Step 2: Configure Visualization Settings")
            
            with st.expander("Pangenome Category Thresholds", expanded=False):
                st.markdown("""
                Adjust the thresholds used to categorize genes in the pangenome:
                - **Core genes**: Present in X% or more of genomes
                - **Soft-core genes**: Present in Y% to X% of genomes
                - **Shell genes**: Present in Z% to Y% of genomes
                - **Cloud genes**: Present in less than Z% of genomes
                """)
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    core_threshold = st.slider("Core gene threshold (%)", 
                                              min_value=90, max_value=100, value=99) / 100
                with col2:
                    softcore_threshold = st.slider("Soft-core gene threshold (%)", 
                                                 min_value=80, max_value=95, value=95) / 100
                with col3:
                    shell_threshold = st.slider("Shell gene threshold (%)", 
                                              min_value=5, max_value=40, value=15) / 100
            
            # Visualization tabs
            st.markdown("### Step 3: View Visualizations")
            tabs = st.tabs([
                "Gene Categories", 
                "Pangenome Distribution", 
                "Gene Frequency", 
                "Rarefaction Curve",
                "Presence/Absence Matrix", 
                "Export Data"
            ])
            
            with tabs[0]:
                st.markdown("#### Gene Category Distribution")
                log_scale = st.checkbox("Use logarithmic scale", key="bar_log")
                plot_bar(results['df'], log_scale)
            
            with tabs[1]:
                st.markdown("#### Pangenome Distribution")
                plot_pie(results['roary'], core_threshold, softcore_threshold, shell_threshold)
            
            with tabs[2]:
                st.markdown("#### Gene Frequency Distribution")
                st.markdown("This plot shows the number of genes present in X genomes.")
                generate_frequency_plot(results['roary'])
            
            with tabs[3]:
                st.markdown("#### Gene Accumulation Curve (Rarefaction)")
                st.markdown("""
                This curve shows how the total number of genes increases as more genomes are added.
                It can help determine if your dataset has sampled the complete pangenome.
                """)
                num_permutations = st.slider(
                    "Number of random permutations", 
                    min_value=10, 
                    max_value=200, 
                    value=50
                )
                generate_rarefaction_curve(results['roary'], num_permutations)
            
            with tabs[4]:
                st.markdown("#### Presence/Absence Matrix")
                
                if 'tree' not in results:
                    st.warning("Tree file (accessory_binary_genes.fa.newick) not provided. Matrix will be displayed without phylogenetic ordering.")
                
                col1, col2 = st.columns(2)
                with col1:
                    cluster_rows = st.checkbox("Cluster rows", value=True)
                with col2:
                    cluster_cols = st.checkbox("Cluster columns", value=False)
                
                plot_matrix(results['roary'], results.get('tree'), cluster_rows, cluster_cols)
            
            with tabs[5]:
                st.markdown("#### Export Processed Data")
                st.markdown("""
                Download the processed data for further analysis in various formats:
                - CSV format (comma-separated)
                - TSV format (tab-separated)
                - Gene frequency statistics
                """)
                export_data(results['roary'])
        
        else:
            if 'df' not in results or results['df'] is None:
                st.error("Failed to process summary statistics file.")
            if 'roary' not in results or results['roary'] is None:
                st.error("Failed to process gene presence/absence matrix.")

    # Add information about the app
    with st.expander("About this app"):
        st.markdown("""
        ### About Roary Output Visualizer
        
        This application is designed to visualize the output from the Roary pangenome analysis pipeline.
        This is created by Tamoghna Das please cite this while using it.
        **Features:**
        - Interactive visualizations of gene distribution
        - Pangenome category analysis
        - Gene frequency plots
        - Rarefaction curve analysis
        - Presence/absence matrix visualization
        - Data export in multiple formats (CSV, TSV, SVG, PNG)
        
        **Notes:**
        - All processing happens in your browser - no data is sent to any server
        - For large datasets, some visualizations may take longer to generate
        - SVG exports are ideal for publication-quality figures
        
        For more information about Roary, visit the [Roary GitHub page](https://github.com/sanger-pathogens/Roary).
        """)

    # Cleanup function
    def cleanup():
        if os.path.exists(st.session_state.temp_dir):
            shutil.rmtree(st.session_state.temp_dir)
    
    # Register cleanup function
    import atexit
    atexit.register(cleanup)

# Run the app
if __name__ == "__main__":
    main()
