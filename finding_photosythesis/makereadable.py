from Bio import Entrez
import os

#basic helper function for accessing gnee information based on an accession code
def fetch_gene_info(identifier):
    Entrez.email = "user@byu.edu"  # Always include your email address
    handle = Entrez.efetch(db="protein", id=identifier, rettype="gb", retmode="text")
    record = handle.read()
    handle.close()
    return record

#use if you want to make large blast files more readable (split up into many different files that open quicker)
def make_readable_len(file_path):
    count = 1
    line_num = 0
    with open(file_path, 'r') as infile:
        for line in infile:
            outpath = "shorter_to_read/" + str(count) + "_" + file_path + ".txt"
            with open(outpath, 'a') as outfile:
                if line_num >= 50000:
                    line_num = 0
                    count += 1
                else:
                    outfile.write(line)
                    line_num += 1

#collect counts of all organisms and hits in a blastn .out file
def collect_organisms(file_path):
    with open(file_path, 'r') as infile:
        outpath = file_path[:-4] + "_organisms.txt"
        with open(outpath, 'w') as outfile:
            orga_count = {}
            switcheroo = False
            for line in infile:
                if switcheroo and line != '\n':
                    if line[0]!=">":
                        l = line.split()
                        if l[1] + " " + l[2] in orga_count:
                            orga_count[l[1]+ " " +l[2]] += 1
                        else:
                            orga_count[l[1]+ " " +l[2]] = 1
                    else:
                        switcheroo = False
                elif len(line.split() )!= 0:
                    if line.split()[0] == "Sequences":
                        switcheroo = True
            for k,v in sorted(orga_count.items(), key=lambda x: x[1]):
                outfile.write(f"{k}: {v}\n")

#Use to collect genes from blastx files
def accession_to_gene_name(file_path):
    with open(file_path, 'r') as infile:
        with open('s1.txt', 'w') as outfile:
            id = ""
            for line in infile:
                if line[0]==">":
                    id = line.split()[0][1:]
                if line[1:6] == "Score":
                    outfile.write(f"{id}, {line.split(',')[0]}, {line.split(',')[1]}"+"\n")
    with open('s1.txt', 'r') as file:
        with open('s2.txt', 'w') as outfile:    
            for identifier in file:
                try:
                    outfile.write(f"${identifier}{fetch_gene_info(identifier.split(',')[0])}")
                except:
                    print(f'{identifier} code access failure')
    with open('s2.txt', 'r') as file:
        with open(file_path[:-4] + "_genelist.txt", 'w') as outfile:
            acce = ""
            for line in file:
                if line[0] == "$":
                    acce = line[1:-1]
                else:
                    l=line.replace(" ", "")
                    if len(l) > 7:
                        if l[:7]=="/gene=\"":
                            outfile.write(f"{acce}, {l[7:-2]}\n")
    os.remove('s1.txt')
    os.remove('s2.txt')

#collect counts of all top scoring significant alignment associated organisms in a blastn search
def collect_winners(file_path):
    with open(file_path, 'r') as infile:
        outpath = 'winners.txt'
        with open(outpath, 'a') as outfile:
            outfile.write(f'Winners of {file_path}'+'\n')
            orga_count = {}
            orga_score = {}
            orga_eval = {}
            switcheroo = False
            for line in infile:
                if switcheroo and line != '\n':
                    if line[0]!=">":
                        l = line.split()
                        if l[1] + " " + l[2] in orga_count:
                            orga_count[l[1]+ " " +l[2]] += 1
                            switcheroo = False
                        else:
                            orga_count[l[1]+ " " +l[2]] = 1
                            switcheroo = False
                        if l[1] + " " + l[2] in orga_score:
                            orga_score[l[1]+ " " +l[2]] += l[-2]
                            switcheroo = False
                        else:
                            orga_score[l[1]+ " " +l[2]] = l[-2]
                            switcheroo = False
                        if l[1] + " " + l[2] in orga_eval:
                            orga_eval[l[1]+ " " +l[2]] += l[-1]
                            switcheroo = False
                        else:
                            orga_score[l[1]+ " " +l[2]] = l[-1]
                            switcheroo = False
                    else:
                        switcheroo = False
                elif len(line.split() )!= 0:
                    if line.split()[0] == "Sequences":
                        switcheroo = True
            for k,v in sorted(orga_count.items(), key=lambda x: x[1]):
                outfile.write(f"{k},{v},{orga_score[k]},{orga_eval[k]}\n")
            outfile.write('\n')

#use after performing multiple of collect winners to get totals of top hits for each section across all files ran
def compile_counts(file_path):
    with open(file_path, 'r') as infile:
        outpath = 'winner.txt'
        with open(outpath, 'a') as outfile:
            outfile.write(f'Winners of {file_path}'+'\n')
            orga_count = {}
            for line in infile:
                try:
                    n_and_c = line.split(': ')
                    if n_and_c[0] in orga_count:
                        orga_count[n_and_c[0]] += int(n_and_c[1])
                    else:
                        orga_count[n_and_c[0]] = int(n_and_c[1])
                except:
                    continue
            for k,v in sorted(orga_count.items(), key=lambda x: x[1]):
                outfile.write(f"{k}: {v}\n")
            outfile.write('\n')

#use to take figure 1 postoutput to usable input for the figure in the form of a csv file
def figure_one_csv(file_path):
    with open(file_path, 'r') as infile:
        with open("fig1.csv", 'w') as outfile:
            gene_val = "psaA"
            outfile.write("gene,alignment,Score,E-value\n")
            for line in infile:
                if line[0:3] == "psa":
                    gene_val = line[0:-1]
                else:
                    l = line.split()
                    outfile.write(f"{gene_val},{l[0]+l[1]},{l[2]},{l[3]}\n")
    
#use to take figure 1 postoutput (from accession to gene name) to usable input for the figure in the form of a csv file
def figure_three_csv(*file_names):
    with open("fig3.csv", 'w') as outfile:
        outfile.write("\"accid\",\"score\",\"expect\",\"gene_name\"\n")
    for file_name in file_names:
        with open(file_name, 'r') as infile:
            with open("fig3.csv", 'a') as outfile:
                for line in infile:
                    l = line.split(", ")
                    if l[3][0:3] == "pet" or l[3][0:3] == "psa" or l[3][0:3] == "psb" or l[3][0:3] == "atp":
                        outfile.write(f"\"{l[0]}\",\"{l[1].split()[2]}\",\"{l[2].split()[2]}\",\"{l[3][:-1]}\"\n")