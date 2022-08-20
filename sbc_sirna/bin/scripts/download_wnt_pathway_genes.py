import requests

# call maayanlab harmonizome API to get wnt pathway info
resp = requests.get("https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/Wnt+signaling+pathway/PANTHER+Pathways", 
		headers={"Content-Type" : "application/json"})

# exract the pathway members
targets = resp.json()['associations']

tnames = []
# tthresh = []
for item in targets:
	tnames.append(item['gene']['symbol'])
	# tthresh.append(item['thresholdValue'])

# write this list out for comparison with correlated genes
with open('../../data/wnt_pathway_genes.txt', 'w') as f:
	f.write('\n'.join(tnames))
	f.write('\n')
