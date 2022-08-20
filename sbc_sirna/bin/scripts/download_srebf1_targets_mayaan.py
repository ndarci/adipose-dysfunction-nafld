import requests

# call maayanlab harmonizome API to get SREBF1 info
resp = requests.get("https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/SREBF1/TRANSFAC+Curated+Transcription+Factor+Targets", 
		headers={"Content-Type" : "application/json"})

# exract the downstream targets
targets = resp.json()['associations']

tnames = []
# tthresh = []
for item in targets:
	tnames.append(item['gene']['symbol'])
	# tthresh.append(item['thresholdValue'])

# write this list out for comparison with DE genes
with open('../../data/srebf1_targetgenes_mayaanlab.txt', 'w') as f:
	f.write('\n'.join(tnames))
	f.write('\n')
