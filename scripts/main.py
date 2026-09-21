#!/usr/bin/python3

from os import listdir as ls
from os import chdir as cd
from os import getcwd as pwd
from core.general import *
import subprocess

cmd = '''awk '{if ($3 > max) {max=$3; ref=$1; t=t + 1}} END {print t * 2 - 1 "," t * 2 "p"}' barcode07_sorted_idxstats.txt | xargs -iR sed -n R v3_reference.fa > ev_ref.fa'''

def read_from_file(name:str) -> str:
	with open(name,'r') as f:
		return ''.join(f.readlines())
	return ''

debug(pwd())
template = read_from_file('core/general.sh')
reference = 'v3_reference.fa'
pre = []

def read_barcode(x: str) -> int:
	return x[:x.find('.fastq')] if '.fastq' in x else x
	x = x[x.find('barcode') + len('barcode'):]
	ans = ''
	for i in x:
		if i == '.':
			break
		ans = ans + i
	return ans

def get_barcode():
	for i in [i for i in ls() if i != 'result' and i != 'final']:
# 		debug(f'Parsing {i}')
		try:
			cd(f'{i}/demulti')
		except:
			continue
		for i in [[read_barcode(j),i,f'{j}'] for j in ls() if read_barcode(j) != '']: pre.append(i);
		cd(root)

def make_script(i, j = None):
	barcode = i[0]
	dir = i[1]
	fastq = i[2]
	now = template
	if not exist(f'{pwd()}/core/lib/references/{barcode}.fa'):
		warning(f'{pwd()}/core/lib/references/{barcode}.fa does not exist, using the general one.')
	now = now.replace('${reference}', f'references/{barcode}.fa' if exist(f'{pwd()}/core/lib/references/{barcode}.fa') else reference)
	now = now.replace('${barcode}',str(barcode) if j is None else f'{barcode}-{j}')
	now = now.replace('${fastq}',fastq)
	now = now.replace('${name}','auto')
	if j is not None:
		now = now.replace('$selected', str(j))
	return now

def fake_cp(i):
	barcode = i[0]
	dir = i[1]
	fastq = i[2]
	selected = i[3]
	cd(root)
	cd(dir)
	os.system(f'mkdir "{root}/result/{barcode}-{selected}"')
	for i in ls():
		os.system(f'ln "{root}/{dir}/{i}" "{root}/result/{barcode}-{selected}/{i}" -s')
	cd(root)
	return

get_barcode()
info('List of operation:')
for i in pre:
	
	ref = ('REF' if exist(f'{pwd()}/core/lib/references/{i[0]}.fa') else 'ORG')
	warning('\t'.join([ref] + [str(j) for j in i]))

if input('Continue operation? [y/N]').lower() != 'y': info('Operation Canceled');quit();

info('Trying to make result directory')
cd(root)
os.system('mkdir result')

for i in pre:
	barcode = i[0]
	dir = i[1]
	fastq = i[2]
	i.append(0)
	for j in range(1, 5 + 1):
		i[3] = j
		cd(root)
		try:
			read_from_file(f'result/{barcode}-{j}/k2.report.txt')
			warning(f'barcode {barcode}-{j} has already exists. skipping.')
			continue
		except:
			pass
		info(f'Making directory for barcode {barcode}-{j} name {fastq}')
		fake_cp(i)
	#	os.system(f'cp -r {dir} result/{barcode}')
		cd(f'result/{barcode}-{j}')
#		info(make_script(i, j))
#		exit()
		subprocess.call('bash -c ". "/home/ysh/miniconda3/etc/profile.d/conda.sh"; conda activate summer; $sz"', shell = True, env = {'sz': make_script(i, j)})

cd(root)

file_to_sort = [['webpage'],['bam'],['coverage'],['species']]

for i in file_to_sort:
	info(f'Making directory {i[0]}')
	os.system(f'mkdir {i[0]}')

for i in pre:
	barcode = i[0]
	dir = i[1]
	fastq = i[2]
	for j in range(1, 5 + 1):
		barcode = i[0] + '-' + f'{j}'
		info(f'Moving webpage for barcode {barcode}')
		os.system(f'cp result/{barcode}/taxonomy.krona.html webpage/{barcode}.html')

		info(f'Moving bam file for barcode {barcode}')
		os.system(f'cp result/{barcode}/barcode{barcode}_sorted.bam bam/{barcode}.bam')
		os.system(f'cp result/{barcode}/barcode{barcode}_sorted.bam.bai bam/{barcode}.bam.bai')

		info(f'Moving coverage file for barcode {barcode}')
		os.system(f'cp result/{barcode}/barcode{barcode}_coverage.txt coverage/{barcode}_coverage.txt')

		info(f'Moving species file for barcode {barcode}')
		os.system(f'cp result/{barcode}/k2.report.txt species/{barcode}.txt')

info('Sorting directory')
os.system('mkdir final')

for i in file_to_sort:
	os.system(f'mv {i[0]} final/')

quit()
