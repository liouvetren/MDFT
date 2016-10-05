# import requests

# url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/12312/record/SDF/?record_type=3d"

# with open('test.sdf','w') as f:
# 	r=requests.get(url)
# 	f.write(r.content)


#---------------------------------------------
# make requests using CSV files
#
#---------------------------------------------

import sys
import os
import requests
import subprocess
import sqlite3
import csv
import urllib

def db_create(db):
	conn = sqlite3.connect(db)
	cur = conn.cursor()
	cur.execute('''
		create table if not exists Chemicals(
			id integer primary key AUTOINCREMENT,
			name text not null,
			cid integer unique not null
		)
		''')
	conn.commit()
	conn.close()

class Found(Exception):
	pass

def find_cid_DB(record,cur,buff):
	for name in record:
		cur.execute('''
			select cid from Chemicals where name=?''',(name,))
		tp = cur.fetchone()
		if tp is not None:
			buff[0] = tp[0]
			raise Found

def find_cid_PUG(record,buff):
	api = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
	for name in record:
		url = api+"%s/cids/txt?name_type=complete" % urllib.quote(name)
		r = requests.get(url)
		if r.status_code is 200:
			buff[0] = int(r.content.splitlines()[0])
			raise Found

def update_DB(record,cid,conn):
	cur = conn.cursor()
	cur.execute('''
		select count(*) from Chemicals where cid=?''',(cid,))
	if cur.fetchone()[0] is 0:
		cur.execute('''
			insert into Chemicals(name,cid)
				values(?,?)''',(record[0],cid))
		conn.commit()

def parse(content, iden, ty=str):
	st = content.index(iden)
	st = content.index('\n',st+1)
	st = content.index('\n',st+1)
	ed = content.index('%',st+1)
	return [ty(cell) for cell in content[st:ed].split()]

def createTop(name,cid,ff):
	path = '%s/%d/' % (ff,cid)
	subprocess.call(['mkdir','-p',path])

	top = path+'lig.top'
	if os.path.isfile(top):
		return

	# sdf
	sdf = path+'lig.sdf'
	if not os.path.isfile(sdf):
		url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound' \
		'/cid/%d/record/SDF?record_type=3d' % cid
		r = requests.get(url)
		if r.status_code is not 200:
			print '---> cannot create sdf file'
			return False
		else:
			with open(sdf,'w') as f:
				f.write(r.content)
	# prm crd
	prm = path+'lig.leap.prm'
	crd = path+'lig.leap.crd'
	if not os.path.isfile(prm) or not os.path.isfile(crd):
		temp = path+'temp/'
		subprocess.call(['mkdir','-p',temp])
		try:
			subprocess.call('''cd %s;
				               /bin/bash -i -c "pytleap --lig=../lig.sdf --lfrc=gaff;
				               mv lig.leap.prm ../;
				               mv lig.leap.crd ../" 
				               '''%(temp)
				,shell=True,stdout=open(os.devnull,'w'), stderr=subprocess.STDOUT)
		except:
			print '---> cannot create prm or crd'
			return False

	# TOP file
	top = path+'mytop.txt'
	if not os.path.isfile(top) or True: # Update TOP whatsoever

		try:
			natom = 0
			coord = []
			atom =[]
			cg =[]
			typ = []
			LJ_INDX=[]
			LJ12 = []
			LJ6 = []
			eps = []
			sig = []

			with open(crd,'r') as f:
				f.readline()
				natom = int(f.readline())
				coord.extend([ float(cell) for cell in f.read().split() ])

			with open(prm,'r') as f:
				content = f.read()
				atom = parse(content,'%FLAG ATOM_NAME')
				type_name = parse(content,'%FLAG AMBER_ATOM_TYPE')
				cg = parse(content,'%FLAG CHARGE',ty=float)
				typ = parse(content,'%FLAG ATOM_TYPE_INDEX',ty=int)
				LJ12 = parse(content,'%FLAG LENNARD_JONES_ACOEF',ty=float)
				LJ6  = parse(content,'%FLAG LENNARD_JONES_BCOEF',ty=float)
				LJ_INDX = parse(content,'%FLAG NONBONDED_PARM_INDEX',ty=int)
				# unit conversion
				cg = [cell/18.2223 for cell in cg]
				for i in range(0,natom):
					type_id = typ[i]-1
					ntype = max(typ) 
					LJ_id = LJ_INDX[ntype*type_id+type_id]-1
					A = LJ12[LJ_id]
					B = LJ6[LJ_id]
					if A < 1e-1:
						eps.append(0.0)
						sig.append(1.06)
					else:
						eps.append(B**2/(4*A))
						sig.append((A/B)**(1.0/6))
			with open(top,'w') as f:
				f.write('%s\n' % name)
				f.write('%d\n' % natom)
				f.write('Atom \t Type \t Charge \t Eps(kcal/mol) \t Sig(Angstrom) \t X \t Y \t Z\n')
				for i in range(0,natom):
					f.write('%s \t %s \t %f \t %f \t %f \t %f \t %f \t %f\n' % 
						(atom[i], type_name[i], cg[i],eps[i],sig[i],coord[3*i],coord[3*i+1],coord[3*i+2]))
		except:
			print '---> cannot create top'
			return False
	return True

#-------------- main --------------------

if len(sys.argv) is not 4:
	print fail(0)
	exit(1)
roster=str(sys.argv[1])
db=str(sys.argv[2])
ff=str(sys.argv[3])
print 'namelist :  ', roster
print 'Database :  ', db
print 'ff_PATH  :  ', ff


# check database exists
db_create(db)
conn = sqlite3.connect(db)
cur = conn.cursor()

#--------------------------------------
csvreader=csv.reader(open(roster,'r'))
indx = 0
for record in csvreader:
	indx = indx + 1
	try: # to Find Cid
		buff = [None]
		find_cid_DB(record,cur,buff)
		find_cid_PUG(record,buff)	
		print "Not Recognized: %s in line %d" % (record[0],indx)
		continue
	except Found:
		cid = buff[0]
		update_DB(record,cid,conn)
		if not createTop(record[0],cid,ff):
			print "@ %s in line %d with CID %d" % (record[0],indx,cid)
			print '---------------------------------'
			continue
conn.close()
		