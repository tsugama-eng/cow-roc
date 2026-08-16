# COW-ROC: Codon weight and ribosome-occupied codon utilization platform

## 📌 Quick Platform Overview
*   **Section 1 (Codon Optimization Engine):** Engineers and synthesizes open reading frames using pre-calculated or custom codon weight models.
*   **Section 2 (Codon Weight Acquisition Pipeline):** Resolves de novo codon weights from a user-uploaded expression matrix and sequence datasets. The expression matrix supports CSV, TSV, or TXT formats (comma- or tab-separated) and must include the gene ID in the first column.
*   **Section 3 (Ribo-seq-weighted Cumulative Codon Frequency Browser):** Tracks absolute and comparative cross-species Ribo-seq-weighted cumulative codon frequencies across a 61-base window centered on the P-site.
*   **🔒 Data Privacy:** This platform is built 100% in HTML/JavaScript. Your input sequences and expression datasets are processed entirely within your local browser memory, ensuring complete data confidentiality with zero server-side data egress.

---

## 🔄 Core Workflow: Custom Codon Weight Derivation & Optimization
To calculate custom codon weights and apply them directly to codon optimization, you can transfer data directly from **Section 2 to Section 1** using the following pipeline steps:

1.  **Upload datasets in Section 2:** Input your multi-FASTA CDS file and expression profile matrix. If the matrix contains multiple experimental conditions, the platform automatically log-transforms and averages the expression values by gene for regression.
2.  **Run linear regression:** Adjust your expression thresholds and preferred log scaling modes, then click `[Run linear regression]`. A simple linear regression will be applied to calculate codon weights and predict expression values.
3.  **Pipeline transmission to Section 1:** The computed codon weights are automatically transmitted and loaded into Section 1 under the **"★ User-uploaded pipeline model"** dropdown menu option.
4.  **Optimize codons:** Select this newly activated custom model option in Section 1, paste your target sequence, and click `[Run]` to generate a codon-optimized sequence based on your custom codon weights.

---

## 📊 Data Interpretation

### 1. Codon weights (Section 1 & 2 input and output)
Codon weights capture whether specific codons act as positive or negative drivers of transcript or protein abundance. High positive codon weights are prioritized during codon optimization. Codon weights are shaped by species-specific cellular environments of mRNA and protein synthesis and turnover.

### 2. Ribo-seq-weighted cumulative codon frequencies (Section 3 browser)
Ribo-seq-weighted cumulative codon frequencies capture the codon frequencies associated with ribosome-stalling sites. The 61-base window is centered on the P-site, where position 0 corresponds to the P-site. The positions from -5 to 5 correspond to the ribosome-occupied sites, which are visualized in the heatmap using distinct color patterns to differentiate them from the surrounding sites.

*   **Absolute density profile:** Displays a heatmap based on absolute frequency data. Unchecking "Hide stop codons" displays stop codon frequencies. In vertebrates and plants, these frequencies are so dominant at the A-site that they obscure the color variations of other codons.
*   **Species contrast (Δ view):** Displays a heatmap based on differential frequency data between two selected species. High absolute delta values ($|\Delta|$) pinpoint codons with significant frequency differences around the ribosome-stalling sites between those species.
*   **Amino acid aggregation:** Displays a heatmap based on frequency data aggregated by synonymous codons. This option allows users to explore the relevance of amino acids to ribosome stalling rather than the effects of individual codons.

---

## **Scripts and programs used for original workflows (in "scripts")**
These files provide the complete (command) pipeline from raw data processing to final analysis.

| Category | File Name | Description |
| :--- | :--- | :--- |
| **Workflows** | `RNA-Seq_script.txt` | Raw RNA-seq reads to TPM values. |
| | `tpm2demand_script.txt` | CDS to codon demand/supply metrics. |
| | `cds2index.txt` | CDS to codon indices and regression models. |
| | `Ribo-seq_script.txt` | Ribo-seq reads to codon stalling profiles. |
| **Programs** | `count2tpm.pl`, `bindfpkm.pl` | RNA-seq processing utilities. |
| | `cds_codon_analysis.pl`, `codon_demand_analysis.pl` | Codon usage and demand calculators. |
| | `stAI_calc.py`, `run_gtAI.py`, `DtAI_calc.py` | Codon index calculators (stAI, gtAI, DtAI). |
| | `codon_exp_regression.py` | Linear regression model generator (Codon Weights). |
| | `predict_tpm.py`, `join_files.py` | Expression prediction and data merging. |
| | `extract_codon_stalling.py`, `extract_ribo_peaks.py` | Ribo-seq P-site and peak extraction. |
| | `codon_read_pileup.py`, `codon_read_pileup2.py` | Ribo-seq frequency calculation (O/E ratios). |
| **Plotting** | `correlation_analysis.R` | Stats and plots for TPM/HL/PA correlations. |
| | `codon_weight_heatmap_UMAP.R` | Interspecies comparison (Heatmap/UMAP). |
| | `Ribo-Seq_analysis_script.R` | Ribo-seq dynamics and stalling visualizations. |

---

## References


*   **[Main publication]**  
    Tsugama D, Kambara K (2026) **Decoding universal principles of codon-mediated regulation of gene expression.** *bioRxiv*  
    [https://doi.org/10.64898/2026.05.27.728126](https://doi.org/10.64898/2026.05.27.728126)

*   **[Data and code repository]**  
    Tsugama D, Kambara K (2026) **Data and Code for: Decoding universal principles of codon-mediated regulation of gene expression.** *Figshare*  
    [https://doi.org/10.6084/m9.figshare.32414811](https://doi.org/10.6084/m9.figshare.32414811)
