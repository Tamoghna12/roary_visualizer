# Roary Visualizer

![Roary Visualizer Banner](https://via.placeholder.com/1200x300?text=Roary+Visualizer)  

A Streamlit-based interactive application to visualize the output of [Roary](https://sanger-pathogens.github.io/Roary/), a tool for pan-genome analysis. This application provides dynamic visualizations of gene categories, pangenome distribution, gene frequency, rarefaction curves, and presence/absence matrices, along with extensive data export options.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [File Structure](#file-structure)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Overview

Roary Visualizer is designed for researchers and bioinformaticians who need to explore and interpret the output from the Roary pipeline. The tool processes Roary output files (such as `summary_statistics.txt`, `gene_presence_absence.csv`, and an optional Newick tree file) and provides interactive visualizations that help users explore:
- **Gene Category Distribution:** Visualizes the count of genes in various categories (core, soft-core, shell, cloud).
- **Pangenome Distribution:** Analyzes the distribution of genes across multiple genomes.
- **Gene Frequency and Rarefaction Analysis:** Displays how gene presence accumulates with the addition of more genomes.
- **Presence/Absence Matrix:** Offers a heatmap representation of gene distribution across genomes, with optional phylogenetic ordering when a Newick tree is provided.

## Features

- **Interactive Visualizations**: Explore multiple interactive plots such as bar graphs, pie charts, and heatmaps.
- **Custom Thresholds**: Dynamically adjust the thresholds for core, soft-core, shell, and cloud genes.
- **Data Export**: Download visualizations as PNG or SVG files, and export processed data in CSV or TSV format.
- **Phylogenetic Tree Integration**: (Optional) Order and display the presence/absence matrix based on a provided Newick tree file.
- **In-Browser Processing**: All processing occurs client-side; your data stays on your machine.
- **User-Friendly Interface**: Built on Streamlit to provide an intuitive and responsive user experience.

## Installation

### Prerequisites

Make sure you have Python 3.8 or later installed. The following Python libraries are required:
- streamlit
- pandas
- numpy
- matplotlib
- seaborn
- plotly
- biopython

### Setup

1. **Clone the Repository**

   ```bash
   git clone https://github.com/yourusername/roary-visualizer.git
   cd roary-visualizer
