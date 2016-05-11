import os



def kallisto(orderNumber,outFolder):
	mapFolder=outFolder + 'mapped/'
	
	
	try:
		os.mkdir(mapFolder)
	except OSError:pass
	
	folders=[]
	readFolder = '/home/david/py/readInput/reid_' + str(orderNumber) 
	
	lanes=os.listdir(readFolder)
	for lane in lanes:
		files=os.listdir(readFolder + '/' + lane)
		samples=[]
		for file in files:
			if file.split('.')[-1]=='gz':
				sample=file.split('_')[0]
				if sample not in samples:
					samples.append(sample)
		
		for sample in samples:
			filesToUse=[]
			for file in files:
				if file.split('_')[0]==sample:
					filesToUse.append(readFolder + '/' + lane + '/' + file)
			output=sample.replace('S','')
		
			command='kallisto quant --single -l 80 -s 20 -i ~/py/database/Danio_rerio/Danio_rerio.GRCz10.cdna.all.idx -o ' + mapFolder + output + ' ' + readFolder
			for file in filesToUse:
				command = command + ' ' + file
			print command
			os.system(command)
			folders.append(mapFolder + output)
	return folders
	

def map(orderNumber,outFolder):
	if outFolder[-1]!='/':
		outFolder=outFolder + '/'
	try:
		os.mkdir(outFolder)
	except OSError:pass
	
	folders=kallisto(orderNumber,outFolder)
	folders.sort()
	first=True
	
	geneinfo=dict()
	
	for folder in folders:
		thisCount = dict()
		for line in open(folder + '/abundance.tsv'):
			if len(line)>2 and line[0] != 't':
				line=line.strip('\n').split('\t')
				gene=line[0]
				if first:
					geneinfo[gene]=line[1] + '\t' + line[2]
				thisCount[gene]=float(line[4])
		vars()[folder.split('/')[-1]]=thisCount
		
		first=False
	IDNames=dict()
	for line in open('CHDM RNA-seq samples - Sheet1.tsv'): #Google sheet downloaded as tsv
		if len(line)>4 and line[0] != '#':
			line=line.split('\t')
			IDNames[line[0]]=line[1]
			
	write=open(outFolder + outFolder.replace('/','') + '.tsv','w')
	write.writelines('gene\tlength\teff_length')
	for folder in folders:
		sampleID=folder.split('/')[-1]
		name=IDNames[sampleID]
		write.writelines('\t' + name)
	for gene in geneinfo.keys():
		write.writelines('\n' + gene + '\t' + geneinfo[gene])
		for folder in folders:
			sampleID=folder.split('/')[-1]
			write.writelines('\t' + str(vars()[sampleID][gene]))
	write.close()
print 'Type: \'map(orderNumber,outFolder)\''	
		
