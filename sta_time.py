import os
import re
import pandas as pd

keys = ['fastqc', 'spades', 'rnammer', 'prodigal', 'kegg /home/data/functional-database/KEGG/kegg-prokaryotes.dmnd']

d = '.'
df_total = pd.DataFrame()
for tax in os.listdir(d):
	log = os.path.join(d, tax, 'output/log')
	if os.path.isfile(log):
		print(log)
		a = {}
		with open(log, 'rt') as h:
		    for i in h:
		        if re.search('INFO ', i):
		            key = re.search('INFO (.+)', i).group(1)
		            tmp = i.strip().split()
		            t = tmp[0]+' '+tmp[1]
		            a['end'] = t
		            if key in keys:
		                a[key] = t

		df = pd.to_datetime(pd.DataFrame.from_dict(a, orient='index')[0]).sort_values()
		total_tmp = df['end'] - df['fastqc']
		df = df.shift(-1) - df
		df['total'] = total_tmp
		df.name = tax
		df_total = pd.concat([df_total, df], axis=1)

df_total.T.to_csv('sta_time.csv')

