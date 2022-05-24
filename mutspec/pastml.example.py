from pastml.acr import pastml_pipeline

# Path to the table containing tip/node annotations, in csv or tab format
data = "~/Downloads/data.txt"

# Path to the tree in newick format
tree = "~/Downloads/Albanian.tree.152tax.tre"

# Columns present in the annotation table,
# for which we want to reconstruct ancestral states
# (for Albanian data we only have one column, but multiple columns are also allowed)
columns = ['Country']

# Path to the output compressed map visualisation
html_compressed = "~/Downloads/Albanian_map.html"

# (Optional) path to the output tree visualisation
html = "~/Downloads/Albanian_tree.html"

pastml_pipeline(data=data, data_sep=',', columns=columns, name_column='Country', tree=tree,
                html_compressed=html_compressed, html=html, verbose=True)