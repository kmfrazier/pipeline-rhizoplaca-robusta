<h1> Guides to Lichen Metagenome Photosynthesis Pathway Research </h1>
<h2> Guide to Using tblastnscript.sh Script </h2>

- Change photosystemI_fasta_file_path to your downloaded path of your downloaded .fasta file
- Change sl20079_comb.fasta_database_path to your created database
    - Instructions for creating custom database: https://www.ncbi.nlm.nih.gov/books/NBK569841/
- Change outfile_path_and_name to your desired .out file name and path

<h2> Guide to Creating Figure 1 </h2>

1. Go to https://www.uniprot.org/uniprotkb?dir=ascend&query=trebouxia+aggregata+photosystem+I&sort=gene order the list by gene name and download all of the "psa-", Photosystem I genes (should be a single file) and should look similar to this:

```fasta
>tr|A0A6B9VRS9|A0A6B9VRS9_9CHLO Photosystem I P700 chlorophyll a apoprotein A1 OS=Trebouxia lynnae OX=1825957 GN=psaA PE=3 SV=1
MTISPPERETKKVKIVVDRNPTETSFERWSKPGHFSRTLAKGPATTTWIWNLHADAHDFD
SHTTDLEEISRKVFSAHFGQLGIIFIWLSGMYFHGARFSNYEAWLSDPTHIKPSAQVVWP
IVGQEILNGDVGGGFQGIQITSGFFQLWRAGGITSELQLYSTAIGGLVMAGAMFFAGWFH
YHKAAPKLEWFQNVESMLNHHLAGLLGLGSLAWAGHQIHVSLPINKLLDAGVDPKEIPLP
HEFLLNRDLMAQLYPSFAKGLSPFFTVNWNEYSDFLTFKGGLNPITGGLWLTDTAHHHLA
IAVLFLVAGHQYRTNWGIGHSIKEILEAHKGPFTGEGHKGLYEILTTSWHAQLAINLALF
GSLSIIVAHHMYSMPPYPYLATDYGTQLSIFTHHLWIGGFCIVGAGAHGAIFMVRDYDPT
NNYNNLLDRVIRHRDAIISHLNWVCIFLGFHSFGLYIHNDTLSALGRPQDMFSDTAIQLQ
PVFAQFIQNTHFLAPQLTAPNALAGTSATWGGDLVAVGGKVAMMPISLGTADFMVHHIHA
FTIHVTVLILLKGVLYARTSRLIPDKANLGFRFPCDGPGRGGTCQVSGWDHVFLGLFWMY
NALSVAIFHFSWKMQSDVWGTVNETGVSHITGGNFAQSANTINGWLRDFMWAQSSQVIQS
YGSALSAYGLIFLGAHFIWAFSLMFLFSGRGYWQELIESIVWAHNKLKVAPAIQPRALSI
TQGRAVGVAHYLLGGIATTWSFFLARIIAVG
>tr|A0A6B9VNL1|A0A6B9VNL1_9CHLO Photosystem I P700 chlorophyll a apoprotein A2 OS=Trebouxia lynnae OX=1825957 GN=psaB PE=3 SV=1
MATKFPKFSQGLAQDPTTRRIWYGIATAHDFESHDGMTEENLYQKIFASHFGQLAIIFLW
TSSLLFHVAWQGNFPQWGRDPLHVRPIAHAIWDPHFGQPAVEAFTRGGASGPVNISTSGV
...
```
2. Create a database using the sl20079_comb.fasta file using the instructions found on https://www.ncbi.nlm.nih.gov/books/NBK569841/
3. Run the photosystem I fasta file against your created database using the tblastn script
4. Collect all of the lines that (note that this was the archaic method used before creating a much better automated script) this can be easily done by searching for "GN=" pasting the gene name(for example: "psaA") then a newline  and proceeding the line starting with "Sequences producing" and collecting subsequent lines until the ">" character, pasting them into a .txt file. It should look something like this:

```BLAST OUTPUT
psaA
D00723:439:HGFFGBCX3:2:1111:5495:88081 1:N:0:TAGCGCTC+TATCCTCT        79.7    2e-17
D00723:439:HGFFGBCX3:1:1110:16899:87827 1:N:0:TAGCGCTC+TATCCTCT       79.7    2e-17
D00723:439:HGFFGBCX3:1:1114:4107:33097 1:N:0:TAGCGCTC+TATCCTCT        79.7    3e-17
D00723:439:HGFFGBCX3:1:1106:10469:84795 2:N:0:TAGCGCTC+TATCCTCT       79.7    3e-17
D00723:439:HGFFGBCX3:1:2116:4551:14873 1:N:0:TAGCGCTC+TATCCTCT        79.7    3e-17
D00723:439:HGFFGBCX3:2:1112:15474:96824 2:N:0:TAGCGCTC+TATCCTCT       78.6    6e-17
D00723:439:HGFFGBCX3:2:1212:11503:13668 2:N:0:TAGCGCTC+TATCCTCT       78.2    9e-17
D00723:439:HGFFGBCX3:1:2202:2399:67459 1:N:0:TAGCGCTC+TATCCTCT        78.2    9e-17
D00723:439:HGFFGBCX3:1:1111:6405:57758 2:N:0:TAGCGCTC+TATCCTCT        78.2    1e-16
D00723:439:HGFFGBCX3:1:1109:3582:91843 2:N:0:TAGCGCTC+TATCCTCT        78.2    1e-16
D00723:439:HGFFGBCX3:1:1113:8914:91351 1:N:0:TAGCGCTC+TATCCTCT        77.8    1e-16
D00723:439:HGFFGBCX3:2:2207:5985:93226 1:N:0:TAGCGCTC+TATCCTCT        77.8    1e-16
D00723:439:HGFFGBCX3:2:1216:20404:38587 1:N:0:TAGCGCTC+TATCCTCT       77.8    1e-16
D00723:439:HGFFGBCX3:1:2214:17415:92466 1:N:0:TAGCGCTC+TATCCTCT       77.8    1e-16
D00723:439:HGFFGBCX3:2:1101:6408:95549 2:N:0:TAGCGCTC+TATCCTCT        77.8    1e-16
D00723:439:HGFFGBCX3:2:2107:12170:47939 1:N:0:TAGCGCTC+TATCCTCT       77.8    1e-16
D00723:439:HGFFGBCX3:2:2216:6553:91313 1:N:0:TAGCGCTC+TATCCTCT        77.0    3e-16
D00723:439:HGFFGBCX3:1:1202:17052:14190 2:N:0:TAGCGCTC+TATCCTCT       77.0    3e-16
D00723:439:HGFFGBCX3:1:2107:12631:92626 1:N:0:TAGCGCTC+TATCCTCT       77.0    3e-16
D00723:439:HGFFGBCX3:2:2103:11746:95720 2:N:0:TAGCGCTC+TATCCTCT       76.6    3e-16
D00723:439:HGFFGBCX3:1:2114:2559:3875 2:N:0:TAGCGCTC+TATCCTCT         76.6    4e-16
psaB
D00723:439:HGFFGBCX3:2:2107:12170:47939 1:N:0:TAGCGCTC+TATCCTCT       77.8    1e-16
D00723:439:HGFFGBCX3:2:2216:6553:91313 1:N:0:TAGCGCTC+TATCCTCT        77.0    3e-16
D00723:439:HGFFGBCX3:1:1202:17052:14190 2:N:0:TAGCGCTC+TATCCTCT       77.0    3e-16
D00723:439:HGFFGBCX3:1:2107:12631:92626 1:N:0:TAGCGCTC+TATCCTCT       77.0    3e-16
D00723:439:HGFFGBCX3:2:2103:11746:95720 2:N:0:TAGCGCTC+TATCCTCT       76.6    3e-16
D00723:439:HGFFGBCX3:1:2114:2559:3875 2:N:0:TAGCGCTC+TATCCTCT         76.6    4e-16
...
```
5. Once this is done, run the text through the the figure_one_csv function and run the figure1.r file to generate the figure

<h2>Guide to Using blastxscript.sh or blastnscript.sh(skip step 6) Script</h2>

1. In line 7 replace the REPLACE_JOB_NAME with your preferred job name
2. In line 8 replace the sample_user@byu.edu with your preferred email
3. In line 2, 13, and 14 replace all instances of group_name with your file sharing group name
4. In line 26 replace ../../home/USER/groups/group_name/file.fasta with the address of the desired .fasta file
5. In line 26 replace ../../home/USER/groups/group_name/file.out with the address of the desired location of the .out results file
6. In line 26 either remove or replace -taxids 000000 depending on your desire to limit results by taxonomy or not at all
7. In line 26 either remove or replace -evalue 1e-5 with the desired E-value threshhold or if you want all of them
8. Run the file via the command line by using the following commands(replace group_name with group sharing folder): 
    - sg group_name "sbatch blastxscript.sh"
    - sg group_name "sbatch blastnscript.sh"

<h2>Guide to Using/Interpreting/Cleaning BLAST .out Files</h2>

- See comments above each function to understand purpose
- biopython is required to install (pip install biopython)
- Change line 6 to reflect your own email address (required for function to work)

<h2>Guide to Creating Figure 3</h2>

1. Run accession_to_gene_name on all of the blastx output (output should be (blastx output name)_genelist.txt)
2. Run figure_three_csv on the (blastx output name)_genelist.txt (output from previous step)
3. Run figure3.r to generate the figure
4. To generate the KEGG diagram of the photosynthesis pathway, go to https://www.kegg.jp/pathway/map00195 and click the "+" on the "Color" sidebar and paste the text from keggcolorocodes.txt and click "Exec" this should color label the gene boxes.
5. Colidation, text, and legends were added manually via Microsoft Word