This README covers sections that explain the project, its features, setup instructions, usage, contributing guidelines, and more, ensuring the repo looks professional and is easy for users to understand.

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
   ```

2. **Create a Virtual Environment**

   ```bash
   python -m venv venv
   ```

3. **Activate the Virtual Environment**

   - On Windows:

     ```bash
     venv\Scripts\activate
     ```

   - On macOS/Linux:

     ```bash
     source venv/bin/activate
     ```

4. **Install the Required Packages**

   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. **Run the Application**

   Launch the Streamlit app using the following command:

   ```bash
   streamlit run app.py
   ```

2. **Upload Roary Output Files**

   In the web interface, you will be prompted to upload the following files:
   - **summary_statistics.txt** (required): Contains gene category statistics.
   - **gene_presence_absence.csv** (required): The gene presence/absence matrix.
   - **accessory_binary_genes.fa.newick** (optional): A Newick tree file for phylogenetic ordering.

3. **Configure Visualization Settings**

   Adjust the thresholds for gene categorization (core, soft-core, shell, cloud) and explore various visualization tabs:
   - **Gene Categories**
   - **Pangenome Distribution**
   - **Gene Frequency**
   - **Rarefaction Curve**
   - **Presence/Absence Matrix**
   - **Export Data**

4. **Export Options**

   Download your visualizations and processed data in various formats (PNG, SVG, CSV, TSV) directly from the interface.

## File Structure

```plaintext
roary-visualizer/
├── app.py                  # Main Streamlit application file
├── requirements.txt        # Python dependencies
├── README.md               # This readme file
└── assets/                 # (Optional) Directory for images or additional assets
```

## Contributing

Contributions are welcome! If you have ideas for improvements, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b feature/my-new-feature
   ```
3. Commit your changes:
   ```bash
   git commit -am 'Add some feature'
   ```
4. Push to the branch:
   ```bash
   git push origin feature/my-new-feature
   ```
5. Open a pull request describing your changes.

For any major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgements

- [Roary](https://sanger-pathogens.github.io/Roary/) for providing a robust pan-genome analysis pipeline.
- The developers and maintainers of Streamlit and the Python libraries used in this project.
- The open-source community for contributions and support.

*Happy visualizing!*

---
