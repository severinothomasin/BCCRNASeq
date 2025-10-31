from pathlib import Path
import scanpy as sc
import scanpy.external as sce
import anndata as ad

def qc_and_filter(adata):
    """
    Annotate mitochondrial, ribosomal, and hemoglobin genes,
    compute QC metrics, and filter low-quality cells/genes.
    """
    # --- Annotate gene categories ---
    var_names_upper = adata.var_names.str.upper()

    # Mitochondrial genes
    adata.var["mt"] = var_names_upper.str.startswith("MT-")

    # Ribosomal genes (both large RPL and small RPS subunits)
    adata.var["ribo"] = var_names_upper.str.startswith(("RPS", "RPL"))

    # Hemoglobin genes (common HB prefixes)
    adata.var["hb"] = var_names_upper.str.startswith(("HB", "HBA", "HBB", "HBD", "HBG"))

    # --- Compute QC metrics ---
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # --- Filter out poor-quality cells and genes ---
    # Basic gene/cell filters (tune as needed)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Filter high mitochondrial percentage
    adata = adata[adata.obs.pct_counts_mt < 10].copy()

    # Optional: filter cells with too few or too many counts
    adata = adata[
        (adata.obs.total_counts > 500) &
        (adata.obs.total_counts < 40000)
    ].copy()

    # Optionally: filter extreme ribosomal or hemoglobin content
    # These thresholds depend on tissue context — you can tune them
    adata = adata[adata.obs.pct_counts_ribo < 60].copy()
    adata = adata[adata.obs.pct_counts_hb < 1].copy()

    return adata  


def run_scrublet_per_sample(adata, sample_key="sample", expected_doublet_rate=0.06):
    # Works on raw counts in .X; ensure not normalized yet
    adata.obs["doublet"] = False
    for s in adata.obs[sample_key].unique():
        print(f'running scrunlet for sample: {s}')
        idx = adata.obs[sample_key] == s
        ad_s = adata[idx].copy()
        sc.pp.scrublet(
            ad_s,
            expected_doublet_rate=expected_doublet_rate,
            n_prin_comps=30,
            sim_doublet_ratio=2.0,
            verbose=False,
        )

        # ad_s.obs["predicted_doublet"] is created by scrublet
        adata.obs.loc[idx, "doublet"] = ad_s.obs["predicted_doublet"].values
        # (Optional) store scores
        adata.obs.loc[idx, "doublet_score"] = ad_s.obs["doublet_score"].values
    return adata


if __name__ == '__main__':

    sc.set_figure_params(dpi=300)

    BASE_DIR = Path("datasets/GSE171892")

    # SAMPLES = [
    #     BASE_DIR / "human_CS12_spinalcord_rep1",
    #     BASE_DIR / "human_CS12_spinalcord_rep2",
    #     BASE_DIR / "human_CS14_brachial_rep1",
    #     BASE_DIR / "human_CS14_thoracic_rep1",
    #     BASE_DIR / "human_CS17_brachial_rep1",
    #     BASE_DIR / "human_CS17_brachial_rep2",
    #     BASE_DIR / "human_CS17_thoracic_rep1",
    #     BASE_DIR / "human_CS17_thoracic_rep2",
    #     BASE_DIR / "human_CS19_brachial_rep1",
    #     BASE_DIR / "human_CS19_brachial_rep2",
    #     BASE_DIR / "human_CS19_thoracic_rep1",
    #     BASE_DIR / "human_CS19_thoracic_rep2"
    # ]

    # datasets = []
    # for sample_path in SAMPLES:
    #     adata = sc.read_10x_mtx(sample_path, var_names="gene_symbols", cache=True)
    #     adata.obs["sample"] = sample_path.name  # Extracts the folder name (e.g. "human_CS12_spinalcord_rep1")
    #     adata.obs_names = [f"{sample_path.name}_{bc}" for bc in adata.obs_names]
    #     datasets.append(adata)

    # adata = ad.concat(
    #     datasets,
    #     join="inner",
    #     label="sample",    
    #     keys=[a.obs["sample"].iloc[0] for a in datasets]
    # )
    
    # adata = qc_and_filter(adata)

    # sc.pl.violin(adata, ['n_genes_by_counts','total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

    # adata = run_scrublet_per_sample(adata, sample_key="sample", expected_doublet_rate=0.06)
    # adata = adata[~adata.obs["doublet"]].copy()

    # adata.layers["counts"] = adata.X.copy()

    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)

    # sc.pp.highly_variable_genes(
    #     adata,
    #     flavor="seurat_v3",
    #     n_top_genes=2000,
    #     batch_key="sample",
    #     layer="counts"    
    # )
    
    # sc.pl.highly_variable_genes(adata)

    # sc.pp.scale(adata, max_value=10)
    # sc.tl.pca(adata, svd_solver='arpack')

    # sce.pp.harmony_integrate(adata, key='sample')

    # sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    # sc.tl.umap(adata)
    # sc.tl.leiden(adata, resolution=2)

    # adata.write(BASE_DIR / "GSE171892_human_integrated.h5ad")

    adata = sc.read(BASE_DIR / "GSE171892_human_integrated.h5ad")

    cell_populations = {
        # === Progenitors ===
        "RP": ["SOX2", "LMX1A", "MSX1", "MSX2", "PAX3", "WNT1"],
        "FP": ["SOX2", "FOXA2", "FERD3L", "ARX", "SHH", "LMX1B", "NKX6-1"],
        "Progenitor": ["SOX2"],
        "Neuron": ["STMN2","MAP2","ELAVL3"],
        "Oligodendrocyte": ["SOX10", "OLIG2"],
        "Blood": ["SOX17"],
        "Hematopoietic": ["FERMT3"],
        "Erythropoietic": ["KLF1"],
        "Erythrocytes": ["HEMGN"],
        "Erythrocytes_II": ["CA2"],

        # === Mesoderm ===
        "Mesoderm_I": ["FOXC1"],
        "Mesoderm_II": ["FOXC2"],
        "Mesoderm_III": ["TWIST1"],
        "Mesoderm_IV": ["TWIST2"],
        "Mesoderm_V": ["MEOX1"],
        "Mesoderm_VI": ["MEOX2"],

        # === Myoblast and crest progenitors ===
        "Myoblast": ["MYOG"],
        "Neural_crest_progenitor": ["SOX10", "SOX2"],
        "DRG_Progenitor": ["SOX10"],
        "Sensory_neuron_progenitor": ["NEUROD1", "NEUROG1", "NEUROG2"],

        # === Skin ===
        "Skin": ["KRT8"]
    }

    sc.pl.dotplot(adata, cell_populations, groupby='leiden', standard_scale='var')

    cluster_dict = {
    '0': '', '1': '', '2': '', '3': '', '4': '', '5': '', '6': '', '7': '', '8': '', '9': '',
    '10': '', '11': '', '12': '', '13': 'BCC', '14': '', '15': '', '16': 'MN', '17': '', '18': '', '19': 'Roofplate',
    '20': '', '21': '', '22': '', '23': 'BCC', '24': '', '25': '', '26': '', '27': 'Floorplate', '28': '', '29': '',
    '30': '', '31': 'Blood', '32': '', '33': '', '34': '', '35': 'Skin', '36': '', '37': '', '38': '', '39': '',
    '40': '', '41': '', '42': '', '43': ''}

    adata.obs['celltype'] = adata.obs['leiden'].map(cluster_dict)

    sc.pl.umap(adata, color=['celltype'])

    sc.pl.umap(adata, color=['leiden'])
