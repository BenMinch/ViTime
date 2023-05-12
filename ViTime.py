####ViTime####

import os, argparse, re, sys, subprocess
#check python version
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='ViTime: a tool for visualizing Virus time-series data')
parser.add_argument('-i', '--input', help='Input Genome Folder', required=True)
parser.add_argument('-o', '--output', help='Output Folder', required=True)
parser.add_argument('-r', '--Reads', help='Folder of reads', required=True)
parser.add_argument('-t', '--threads', help='Number of threads', required=True)
parser.add_argument('-minid', '--minid', help='Mapping minid 0-100 which is ANI', required=True)
parser.add_argument('-tax', '--taxonomy', help='A tab separated file with at least a column for query(genome names without extension) and order', required=True)
parser.add_argument('-ref', '--reference', help='A csv file with read name corresponding to date (Must have columns SRR_run and Date as well as Sample_Name) but can also have environmental variables', required=True)
parser.add_argument('-envs', '--envs', help='A comma separated list of all the environmental variables you want to test', required=False)
args = parser.parse_args()

input_folder = args.input
output_folder = args.output
reads_folder = args.Reads
threads = args.threads
minid = args.minid
taxonomy = args.taxonomy
reference = args.reference
envs = args.envs
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

###Step 1: Map reads to genomes
print('Mapping reads to genomes')
os.mkdir(output_folder+'/Mapping')
coverm='python scripts/CoverM_Genominator_minimap.py -i '+input_folder+' -r '+reads_folder+' -t '+threads+' -minid '+minid+' -o '+output_folder+'/Mapping'
subprocess.call(coverm, shell=True)

###Step 2: Make a matrix with the average counts for each day
print('Creating read matrix')
os.mkdir(output_folder+'/Matrix')
os.mkdir(output_folder+'/Matrix/renamed')
os.mkdir(output_folder+'/Matrix/combined')
gemap='python scripts/GeMAP.py -c '+output_folder+'/Mapping -r '+output_folder+'/Matrix/renamed -ref '+reference+' -com '+output_folder+'/Matrix/combined -o '+output_folder+ ' -col 3'
subprocess.call(gemap, shell=True)


move= 'mv '+ output_folder+'.csv '+output_folder+'/Matrix'
subprocess.call(move, shell=True)
remove= 'rm -r '+output_folder+'/Matrix/renamed'
subprocess.call(remove, shell=True)
remove2= 'rm -r '+output_folder+'/Matrix/combined'
subprocess.call(remove2, shell=True)

os.mkdir(output_folder+'/Matrix/renamed')
os.mkdir(output_folder+'/Matrix/combined')
gemap2='python scripts/GeMAP.py -c '+output_folder+'/Mapping -r '+output_folder+'/Matrix/renamed -ref '+reference+' -com '+output_folder+'/Matrix/combined -o coverage -col 1'
subprocess.call(gemap2, shell=True)

move= 'mv coverage.csv '+output_folder+'/Matrix'
subprocess.call(move, shell=True)
remove= 'rm -r '+output_folder+'/Matrix/renamed'
subprocess.call(remove, shell=True)
remove2= 'rm -r '+output_folder+'/Matrix/combined'
subprocess.call(remove2, shell=True)

os.mkdir(output_folder+'/Matrix/renamed')
os.mkdir(output_folder+'/Matrix/combined')
gemap2='python scripts/GeMAP.py -c '+output_folder+'/Mapping -r '+output_folder+'/Matrix/renamed -ref '+reference+' -com '+output_folder+'/Matrix/combined -o length -col 4'
subprocess.call(gemap2, shell=True)

move= 'mv length.csv '+output_folder+'/Matrix'
subprocess.call(move, shell=True)
remove= 'rm -r '+output_folder+'/Matrix/renamed'
subprocess.call(remove, shell=True)
remove2= 'rm -r '+output_folder+'/Matrix/combined'
subprocess.call(remove2, shell=True)

###Step 3: Categorize viruses
print('Categorizing viruses')
os.mkdir(output_folder+'/Categorize')
coverage_matrix = pd.read_csv(output_folder+'/Matrix/coverage.csv')
lenth=pd.read_csv(output_folder+'/Matrix/length.csv')
#set length column equal to the second column of length.csv
coverage_matrix['Length']=lenth.iloc[:,1]

#transform all the values by dividing by the length of the genome in that row
#remove unmapped row
coverage_matrix=coverage_matrix[coverage_matrix['Genome']!='unmapped']
#Remove any column that has the word "Unnamed" in it
coverage_matrix=coverage_matrix.loc[:, ~coverage_matrix.columns.str.contains('^Unnamed')]

coverage_matrix.iloc[:,1:]=coverage_matrix.iloc[:,1:].div(coverage_matrix['Length'], axis=0)

#drop the length column
coverage_matrix=coverage_matrix.drop(columns=['Length'])
#Count the total number of columns
total_columns=len(coverage_matrix.columns)
#Create a column that counts how many columns in a row have a value greater than 0.25
coverage_matrix['Count']=coverage_matrix.iloc[:,1:].gt(0.25).sum(axis=1)

#Create a column called Category
coverage_matrix['Category']=''
#loop throught the coverage matrix and categorize each virus
for i in range(len(coverage_matrix)):
    if coverage_matrix['Count'][i]>=0.85*total_columns:
        coverage_matrix['Category'][i]='Persistent'
    if coverage_matrix['Count'][i]< 0.85*total_columns and coverage_matrix['Count'][i]> 0.15*total_columns:
        coverage_matrix['Category'][i]='Occasional'
    if coverage_matrix['Count'][i]<= 0.15*total_columns:
        coverage_matrix['Category'][i]='Sporadic'

#make a bar graph of the counts of each category and save it to the output folder but make x axis labels smaller and horizontal
plt.figure(figsize=(10,10))
plt.xticks(rotation=90, fontsize=8)
coverage_matrix['Category'].value_counts().plot(kind='bar')
plt.savefig(output_folder+'/Categorize/Category_Counts.png')
###Step 4: Calculate Celebrity Score
print('Calculating Celebrity Score')

rpkm_matrix= pd.read_csv(output_folder+'/Matrix/'+output_folder+'.csv')
df=rpkm_matrix
#drop the first column
df2=df.drop('Genome', axis=1)
#Celebrity score
df2['Celeb_score']=0

def super_celeb(df2):
    for i in range(len(df2)):
    #calculate mean for the row
        celeb_score=0
        mean= df2.iloc[i,].mean()
        std_error= df2.iloc[i,:len(df2.columns)-1].std()/np.sqrt(len(df2.columns))
        CI= 1.96*std_error
    #iterate through each column
        for j in range(len(df2.columns)-1):
        #check if value is greater than 1
            if df2.iloc[i,j] > 1:
                if df2.iloc[i,j] > df2.iloc[i,j+1]+CI+mean:
                    celeb_score+=1
        df2['Celeb_score'][i]= celeb_score
    return df2
super_celeb(df2)
df2.to_csv('celeb.csv', index=False)

#send it to celebrity folder
os.mkdir(output_folder+'/Celebrity')
rscript= 'Rscript scripts/celebrity.r '+output_folder+'/Celebrity celeb.csv'
subprocess.call(rscript, shell=True)

##Plots with R script
print('Creating plots')
taxonomy= pd.read_csv(taxonomy, sep='\t')
rpkm_matrix_R= rpkm_matrix
rpkm_matrix_R['Order']=rpkm_matrix_R['Genome'].map(taxonomy.set_index('query')['order'])

rpkm_matrix_R.to_csv('rpkm_matrix_R.csv', index=False)
#Environmental plots
if envs: 
    environment_data= pd.read_csv(reference)

#only select columns that are in the envs list and date column
    environment_data=environment_data[['Date']+envs.split(',')]
#group by date
    environment_data=environment_data.groupby('Date').mean()
#reset index
    environment_data=environment_data.reset_index()
#write to csv
    environment_data.to_csv('environment_data.csv',index=False)

    subprocess.call('Rscript scripts/vitime_plots.r rpkm_matrix_R.csv '+output_folder+' environment_data.csv', shell=True)
else:
    subprocess.call('Rscript scripts/vitime_plots_alt.r rpkm_matrix_R.csv '+output_folder, shell=True)


move= 'mv environment_data.csv '+output_folder
subprocess.call(move, shell=True)
move2= 'mv rpkm_matrix_R.csv '+output_folder
subprocess.call(move2, shell=True)
move3= 'mv celeb.csv '+output_folder
subprocess.call(move3, shell=True)
move4= 'mv Celebrity_membership.png '+output_folder
subprocess.call(move4, shell=True)
