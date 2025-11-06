import pandas as pd
import requests

taxon_id = "9606"  # Human
# go_term = '0009986' # Cell surface
# go_term = '0048870' # Cell motility
# go_terms = ['0048870', '0005925', '0097581', '0030175']
# go_terms=['0003700'] # transcription factors
# go_terms=['0030217'] # T cell
# url_genes_by_go = f'https://api.geneontology.org/api/bioentity/function/GO%3A{go_term}/genes?taxon=NCBITaxon%3A{taxon_id}&relationship_type=involved_in&rows=999999'
# print(url_genes_by_go)
# go_terms=['0016020'] # cell membrane
# go_terms=['0005576'] # ECM
# go_terms=['0003700'] # transcription factors
# go_terms=['0005737'] # cytoplasm
# go_terms=['0005840'] # ribosome
# go_terms=['0005739'] # mitochondria
# negative regulation of gene expression GO:0031047

# specify the go id you want here
go_terms = ["0016020", "0003700", "0005737"]
descriptions = ["cell membrane", "TF", "cytoplasm"]

df_out = None
for i, gt in enumerate(go_terms):
    try:
        url_genes_by_go = f"https://api.geneontology.org/api/bioentity/function/GO%3A{gt}/genes?taxon=NCBITaxon%3A{taxon_id}&relationship_type=involved_in&rows=999999"
        response = requests.get(url_genes_by_go).json()["associations"]
        print(f"Got {len(response)} entries.")
        print(response[0])
        df_out_next = pd.json_normalize(response)
        df_out_next["go_name"] = descriptions[i]
        if df_out is None:
            df_out = df_out_next
        else:
            df_out = pd.concat([df_out, df_out_next], axis=0)
    except requests.exceptions.HTTPError as error:
        print(error)

df_out.to_csv("gene_list.csv")

