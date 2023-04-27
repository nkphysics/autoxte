import subprocess as sp
import os
import pandas as pd
import time
from astroquery.heasarc import Heasarc
from astropy.table import Table
from astropy.visualization import time_support
from astropy.visualization import quantity_support
from astropy.time import Time
import datetime

st = True
observations = []
ras = []
decs = []
cNums = []
prnbs = []

print('##############################################################################################')
print('#                                   AutoXTE Data Reducer V4.1                                #')
print('##############################################################################################')
obj = str(input('Target: '))

q_set = str(input('Write Output Que: [n] '))
q_path = 0
q_name = 0
if q_set == 'y':
	q_path = str(input('Input Que: '))
	if q_path[0] == r"'" or que_path[0] == r'"':
		q_path = q_path.replace("'", "")
		q_path = q_path.replace(" ", "")
		q_path = q_path.replace('"', '')
	else:
		pass
else:
	pass

heasarc = Heasarc()
quantity_support()
time_support()
pca = heasarc.query_object(obj, mission='xtemaster', resultmax=1000000) # calls XTE master catalogue for crab
pca = Table(pca).to_pandas()

cnt = 0
for i in pca['OBSID']: # converts the form of the NICER obsid's to strings
	i = i.decode()
	pca.loc[cnt, 'OBSID'] = str(i)
	cnt = cnt + 1
cnt = 0	
for i in pca['TIME']: # converts times from mjd to datetime format
	t0 = Time(i, format='mjd').to_datetime()
	pca.loc[cnt, 'TIME'] = t0
	cnt = cnt + 1
pca = pca.loc[pca['EXPOSURE'] != 0]

downCommand = ('wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R '+'\'index*\''+' -erobots=off --retr-symlinks https://heasarc.gsfc.nasa.gov/FTP/xte/data/archive/AO')
# beginning of obsid inputs
while st == True: 
	enter = str(input('autoXTE >> '))
	if enter == 'done' or enter == 'Done':
		st = False
	elif enter == 'version' or enter =='-v':
		print('autoXTE V4.1 Date Modified September 14, 2021')
	elif enter == 'sel' or enter == 'Sel' or enter == 'SEL':
		print('Observations Selected:')
		for i in observations:
			print(i)
	elif enter == 'back' or enter == 'Back' or enter == 'BACK':
		print('Old Que: ' + str(observations))
		del observations[-1]
		del ras[-1]
		del decs[-1]
		del cNums[-1]
		del prnbs[-1]
		print('New Que: ' + str(observations))
	else:
		if len(str(enter)) == 5:
			rows = pca.loc[pca['PRNB'] == int(enter)]
			print(rows)
			for i in rows['OBSID']:
				observations.append(i)
			for i in rows['RA']:
				ras.append(i)
			for i in rows['DEC']:
				decs.append(i)
			for i in rows['PRNB']:
				prnbs.append(i)
				prnb = str(i)
				if int(prnb[0]) == 9:
					cycles = [9, 10, 11, 12, 13, 14, 15, 16]
					cNums.append(cycles[int(prnb[1])])
				else:
					cNums.append(int(prnb[0]))
		elif enter == 'all':
			for i in pca['OBSID']:
				observations.append(i)
			for i in pca['RA']:
				ras.append(i)
			for i in pca['DEC']:
				decs.append(i)
			for i in pca['PRNB']:
				prnbs.append(i)
				prnb = str(i)
				if int(prnb[0]) == 9:
					cycles = [9, 10, 11, 12, 13, 14, 15, 16]
					cNums.append(cycles[int(prnb[1])])
				else:
					cNums.append(int(prnb[0]))
		else:
			observations.append(enter)
			row = pca.loc[pca['OBSID'] == enter]
			dt = row['TIME']
			row.reset_index(drop=True, inplace=True)
			dt.reset_index(drop=True, inplace=True)
			# basic if else statement to fix single digit months not having a zero out front
			ras.append(row['RA'][0])
			decs.append(row['DEC'][0])
			prnbs.append(row['PRNB'][0])
			prnb = str(row['PRNB'][0])
			if int(prnb[0]) == 9:
				cycles = [9, 10, 11, 12, 13, 14, 15, 16]
				cNums.append(cycles[int(prnb[1])])
			else:
				cNums.append(int(prnb[0]))
			ras.append(row['RA'][0])
			decs.append(row['DEC'][0])
count = 0
for obsid in observations:
	print('')
	print('--------------------------------------------------------------')
	print('             Prosessing OBSID: ' + str(obsid)) 
	print('--------------------------------------------------------------')
	print('Downloading 00 data...')
	sp.call(str(downCommand) + str(cNums[count]) + '//P' + str(prnbs[count]) + '/' + str(observations[count]) + '/ --show-progress --progress=bar:force', shell=True)
	
	base_dir = os.getcwd()
	orb_file = sp.run('ls P' + str(prnbs[count]) + '/' + str(observations[count]) + '/orbit/FP*', shell=True, capture_output=True, encoding='utf-8')
	orb_file = orb_file.stdout.split('\n')
	orb_file.remove('')
	orb_file = str(base_dir) + '/' + str(orb_file[0])
	os.chdir('P' + str(prnbs[count]) + '/' + str(observations[count]) + '/pca')
	sp.call('gzip -d *evt.gz', shell=True)
	evt_files = sp.run('ls *evt', shell=True, capture_output=True, encoding='utf-8')
	evt_files = evt_files.stdout.split('\n')
	evt_files.remove('')
	ccn = 1
	for i in evt_files:
		sp.call('barycorr infile=' + str(i) + ' outfile=bcSE_' + str(obsid) + '_' + str(ccn) + '.evt orbitfiles=' + str(orb_file) + ' refframe=ICRS ra=' + str(ras[count]) + ' dec=' + str(decs[count]), shell=True)
		if q_set == 'y':
			read_q = pd.read_csv(q_path)
			newline = pd.Series(data=[str(base_dir) + '/' + 'P' + str(prnbs[count]) + '/' + str(observations[count]) + '/pca/' + 'bcSE_' + str(obsid) + '_' + str(ccn) + '.evt', 'XTE' + str(obsid)], index=['Input', 'Name'])
			read_q = read_q.append(newline, ignore_index=True)
			read_q.to_csv(q_path, index=False)
		else:
			pass
		ccn = ccn + 1		
	os.chdir(base_dir)
	count = count + 1

print('################################################################################################')
print('Process Done')
