import os,sys,argparse,subprocess,re

#define flags
parser = argparse.ArgumentParser(description='A Script to run coverM genome on a directory of genomes with a lot of reads')
parser.add_argument('-i', '--input', help='Input directory of genomes', required=True)
parser.add_argument('-o', '--output', help='Output directory', required=True)
parser.add_argument('-r', '--reads', help='Reads directory', required=True)
parser.add_argument('-minid', help='Minimum percent identity 0-100', required=True)
parser.add_argument('-t', '--threads', help='Number of threads', required=True)

args = parser.parse_args()

#set variables
input_dir = args.input
output_dir = args.output
reads_dir = args.reads
minID = args.minid
threads = args.threads

#make output directory
for i in os.listdir(reads_dir):
	forward= os.path.join(reads_dir,i)
	reverse= re.sub('_1','_2',forward)
	if forward==reverse:
		print('skip')
	else:
		outpath = os.path.join(output_dir,i)
		outpath2= re.sub('.gz','.coverm',outpath)
		cmd= 'coverm genome -1 ' + forward + ' -2 ' + reverse + ' -d ' + input_dir + ' -p minimap2-sr --min-read-percent-identity '+ minID+' -m covered_bases count rpkm length -t '+threads+ ' --min-covered-fraction 0 > ' + outpath2
		subprocess.call(cmd, shell=True)
          
