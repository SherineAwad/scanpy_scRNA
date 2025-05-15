import scanpy as sc
import sys
import importlib_metadata
import argparse

sys.modules['importlib.metadata'] = importlib_metadata

def read_samples(file_path):
    samples = {}
    with open(file_path, 'r') as file:
        for line in file:
            sample_id = line.strip()  # Get the sample name (e.g., x, y, z)
            if sample_id:  # Ensure there's a valid sample name
                # Construct the filename by appending '_filtered_feature_bc_matrix.h5' to the sample name
                filename = f"{sample_id}_filtered_feature_bc_matrix.h5"
                samples[sample_id] = filename
    return samples



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile')
    parser.add_argument('obj_name') 
    args = parser.parse_args()
    
    inputfile = args.inputfile
    obj_name = args.obj_name 

    adatas = {}
    samples = read_samples(inputfile) 
    print(samples)
    for sample_id, filename in samples.items():
        print(f"Reading {filename}...")  # Optional: Print status
        adata = sc.read_10x_h5(filename)  # Read the file
        adata.var_names_make_unique()
        adata.obs['sample'] = sample_id
        adatas[sample_id] = adata

    combined_adata = sc.concat(adatas.values(), label='sample', keys=adatas.keys())
    combined_adata.var["mt"] = combined_adata.var_names.str.startswith("mt-")

    sc.pp.calculate_qc_metrics(
         combined_adata, qc_vars=["mt"], inplace=True, log1p=True
    )

    sc.pl.violin(
         combined_adata,
         ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
         jitter=0.4,
         multi_panel=True, save="_QC.png"
    )


    combined_adata = combined_adata[
      (combined_adata.obs['n_genes_by_counts'] > 800) &
      (combined_adata.obs['n_genes_by_counts'] < 6000) &
      (combined_adata.obs['total_counts'] > 1200) &
      (combined_adata.obs['total_counts'] < 30000) &
      (combined_adata.obs['pct_counts_mt'] < 25), :
      ]



    sc.pp.filter_cells(combined_adata, min_genes=100)
    sc.pp.filter_genes(combined_adata, min_cells=3)


    sc.pl.violin(
        combined_adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True, save="_AfterQC.png"
    )

    filename = obj_name + ".h5ad"
    combined_adata.write(filename)
if __name__ == "__main__":
    main()
