import pandas as pd 

peak_df = pd.read_csv('peaks_per_gene.csv')
peak_df['GeneID'] = peak_df.Gene.apply(lambda x : ".".join(x.split("+")[0].split(".")[:-1]))
tran_df = pd.read_csv('transcripts_per_gene.csv')
tran_df['GeneID'] = tran_df.Gene

merge_df = pd.merge(peak_df, tran_df, on='GeneID', how='inner')
merge_df.to_csv('peak_transcript_density_comp.csv', index=False)
