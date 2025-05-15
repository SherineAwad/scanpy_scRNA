import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse

sys.modules['importlib.metadata'] = importlib_metadata

parser = argparse.ArgumentParser()
parser.add_argument('myObject')
parser.add_argument('name')

args = parser.parse_args()
myObject =  args.myObject
fname = args.name

newObject = "clustered_" + myObject 

combined_adata = sc.read(myObject)

figure_name = fname + "_clusters.png"
sc.tl.leiden(combined_adata, n_iterations=2)
sc.pl.umap(combined_adata, color=["leiden"], save=figure_name, legend_loc="on data") 

marker_genes  = {
    "MG": ["rlbp1a","rlbp1b","gfap","apoea","apoeb","her6","notch1a","notch1b","aqp4","pax6a","prdx6","slc1a3a","slc1a3b","vim"],
    "Rod": ["insm1a","nr2e3","rho","nrl","crx","guca1b","rom1a","rom1b"],
    "Cones": ["prdm1a","opn1mw1","opn1mw2","opn1mw3","opn1mw4","opn1sw1","opn1sw2","arr3a","thrb","gnat2"], 
    "BC": ["vsx1", "sebox","bhlhe23","cabp5a","cabp5b","vsx2","pcp4a","isl1"] ,
    "AC": ["gad1a","gad1b","gad2","slc6a9","tfap2b","prox1a","pax6a","calb2a","calb2b","pcp4a","elavl3","isl1","chata", "th", "stat3"], 
    "HC": ["ompa","lhx1a","cbln4","calb1","nefla","neflb","nefma","nefmb"], 
    "RGC": ["nefla","neflb","nefma","nefmb","sncga","sncgb","thy1","ebf3a","rbfox3a","rbfox3b","isl1","isl2a","isl2b","pou4f1","pou4f2","pou4f3","rbpms"],  
    "Microglia": ["ptprc","csf2rb","mpeg1.1"], 
    "Progenitors": ["her4.2", "her4.3", "her4.4","dla", "ccnd1","ascl1a","neurod1","sox2","cdh2","atoh7"], 
    "Endothelial": ["tie1"],
    "Pericytes": ["kcnj8", "acta2"],
    "Proliferating cells": ["cdk1"], 
    "RPE": ["rpe65a"],
    "Olignocytes ": ["mbpa"],
    "RPC": ["her12", "fgf19","sfrp2", "gli1"], 
    "Astrocytes":["pax2a", "pax2b","igf2a", "igf2b"]  
    }
figure_name = fname + "markerGenes.png"
# Filter out genes not found in the dataset
filtered_marker_genes = {
    cell_type: [gene for gene in genes if gene in combined_adata.var_names]
    for cell_type, genes in marker_genes.items()
}

# Remove cell types with no valid genes
filtered_marker_genes = {
    cell_type: genes for cell_type, genes in filtered_marker_genes.items() if genes
}

# Save figure
figure_name = fname + "_markerGenes.png"
sc.pl.dotplot(
    combined_adata,
    filtered_marker_genes,
    groupby="leiden",
    standard_scale="var",
    save=figure_name
)


figure_name = fname + "all_marker_genes_feature_plots.png"
# Assuming you have already loaded the data and performed UMAP
# Flatten your marker genes dictionary into a list of all genes
all_genes = []
for cluster, genes in marker_genes.items():
    all_genes.extend(genes)

# Set up a grid of subplots: here we decide the grid size based on the number of genes
n_genes = len(all_genes)
ncols = 3  # You can change this based on how you want to organize the plots (e.g., 3 columns)
nrows = (n_genes + ncols - 1) // ncols  # This calculates the number of rows needed to fit all genes

# Create the figure with the appropriate number of subplots
fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 4))  # Adjust figsize as needed
axes = axes.flatten()  # Flatten the axes array to easily loop over them


# Plot each gene in a separate subplot
for i, gene in enumerate(all_genes):
    if gene in combined_adata.var_names:
        # Plot the gene expression as a feature plot in the appropriate subplot
        sc.pl.scatter(combined_adata, color=gene, title=gene, basis='umap',legend_loc='right margin',  ax=axes[i], show=False)
        axes[i].set_title(f'{gene}')  # Optionally, customize the title

plt.tight_layout()
plt.savefig(figure_name, dbi="300")




#For printing a file for each gene 
for cluster, genes in marker_genes.items():
    for gene in genes:
        if gene in combined_adata.var_names:  # Check if the gene is present in the data
            # Plot using sc.pl.scatter with UMAP embedding
            sc.pl.scatter(combined_adata, color=gene, title=f'{cluster} - {gene}', basis='umap', save=f'{fname}_{gene}.png')

combined_adata.obs_names_make_unique()
combined_adata.write(newObject,compression="gzip")
