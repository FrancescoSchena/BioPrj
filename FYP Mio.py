from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from os import walk
import os

# Looking for the files in the /genomes folder
mypath = os.path.abspath(os.getcwd()) + "/filtered_genomes"
ext = "embl"
files = []
for (dirpath, dirnames, filenames) in walk(mypath):
    for file in filenames:
        if file.split(".")[1] == ext:  # Only files with the right extension
            files.append(file)
    break
if not files:  # if files is empty
    raise Exception("Warning there are no files, check the path!")
    quit()  # no need to proceed

FileNames = []
IntergenomicGC = []
listadicodonitotale = []

files_tot = len(files)
for index, name in enumerate(files):
    print("\n")
    print(f"Analyzing file {index+1} of {files_tot}: {name}")
    print("\n")
    record = SeqIO.read(f'filtered_genomes/{name}', ext)
    sequence = record.seq
    FileNames.append(str(name))
    
    coding_sequences = []  # It will contain the non pseudogenes (reversed if needed)
    intergenomic_sequence = []  # It will contain the sequences out of the locations bounds
    intergenomic_sequence_list = []  # It will contain the locations bounds
    intergenomic_sequence_string = ""  # intergenomic sequence
    
    start = 0

    for feature in record.features:

        locations= feature.location
        s = int(locations.nofuzzy_start)
        e = int(locations.nofuzzy_end)
        st = int(locations.strand)
        
        gene_start = int(locations.nofuzzy_start)

        if gene_start > start:
            ig_sequence = sequence[start:gene_start]
            intergenomic_sequence.append(str(ig_sequence))
            intergenomic_sequence_list.append((start, gene_start))

        new_start = int(locations.nofuzzy_end)

        if new_start != len(sequence):
            start = new_start

        genes = str(sequence[s:e])
        #Reversing complementary strands
        if st == -1:
            sequence_negativestrand= genes[::-1]
            sequence_negativestrand= sequence_negativestrand.replace("A","b")
            sequence_negativestrand= sequence_negativestrand.replace("T","A")
            sequence_negativestrand= sequence_negativestrand.replace("b","T")
            sequence_negativestrand= sequence_negativestrand.replace("C","d")
            sequence_negativestrand= sequence_negativestrand.replace("G","C")
            sequence_negativestrand= sequence_negativestrand.replace("d","G")

            genes = str(sequence_negativestrand)

        #Discarding the pseudogenes
        ok = 0
        if len(genes)%3 != 0:
            ok = 1
        elif genes[:3] != "ATG" and genes[:3] != "GTG" and genes[:3] != "TTG":
            ok = 2
        elif genes[-3:] != "TAA" and genes[-3:] != "TAG" and genes[-3:] != "TGA":
            ok = 3
        elif len(genes) == len(sequence):
            ok = 4
        check5 = False
        if check5:
            for i in range(0, len(genes)-3, 3): # Splitting the sequence in codons discarding the last one
                codon = genes[i: i + 3]
                if codon in ["TAA", "TAG", "TGA"]:
                    ok = 5  # Raises an error if a codon in the middle of the gene is a stop codon

        if ok == 0:
            coding_sequences.append(genes)

        #checking for errors
        check = False
        if check and ok != 0:
            print(locations)
            print(ok)
            print(genes[:3])
            print(genes[-3:])
        if check5 and ok == 5:
            print("Stop codon error:")
            print(f"CDS:\n{genes}")
            
    if start != len(sequence):
        ig_sequence = sequence[start:len(sequence)]
        intergenomic_sequence.append(str(ig_sequence))
        intergenomic_sequence_list.append((start, gene_start))

    intergenomic_sequence_string = "".join(intergenomic_sequence)

    # Intergenomic GC content
    IgGC = intergenomic_sequence_string.count("G") + intergenomic_sequence_string.count("C")
    IgGC_content = (IgGC / len(intergenomic_sequence_string)) * 100
    IntergenomicGC.append(IgGC_content)
    
    coding_sequencewhole = "".join(coding_sequences)
    aminoacids_sequence = Seq(coding_sequencewhole).translate()
    
    
    #Codon counting
    dizionariodicose = {"AAA" : 0,
                        "ATA" : 0,
                        "ACA" : 0,
                        "AGA" : 0,
                        "TTA" : 0,
                        "TCA" : 0,
                        "CAA" : 0,
                        "CTA" : 0,
                        "CCA" : 0,
                        "CGA" : 0,
                        "GAA" : 0,
                        "GTA" : 0,
                        "GCA" : 0,
                        "GGA" : 0,
                        "AAT" : 0,
                        "ATT" : 0,
                        "ACT" : 0,
                        "AGT" : 0,
                        "TAT" : 0,
                        "TTT" : 0,
                        "TCT" : 0,
                        "TGT" : 0,
                        "CAT" : 0,
                        "CTT" : 0,
                        "CCT" : 0,
                        "CGT" : 0,
                        "GAT" : 0,
                        "GTT" : 0,
                        "GCT" : 0,
                        "GGT" : 0,
                        "AAC" : 0,
                        "ATC" : 0,
                        "ACC" : 0,
                        "AGC" : 0,
                        "TAC" : 0,
                        "TTC" : 0,
                        "TCC" : 0,
                        "TGC" : 0,
                        "CAC" : 0,
                        "CTC" : 0,
                        "CCC" : 0,
                        "CGC" : 0,
                        "GAC" : 0,
                        "GTC" : 0,
                        "GCC" : 0,
                        "GGC" : 0,
                        "AAG" : 0,
                        "ACG" : 0,
                        "AGG" : 0,
                        "TTG" : 0,
                        "TCG" : 0,
                        "TGG" : 0,
                        "CAG" : 0,
                        "CTG" : 0,
                        "CCG" : 0,
                        "CGG" : 0,
                        "GAG" : 0,
                        "GTG" : 0,
                        "GCG" : 0,
                        "GGG" : 0}
    
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            for cosa in dizionariodicose:
                if codons == cosa:
                    dizionariodicose[cosa] +=1

    
    aminoacids_sequence = Seq(coding_sequencewhole[3:]).translate()
    boh = len(aminoacids_sequence)

    listadicodoniperfile = [dizionariodicose[i] / boh *100 for i in dizionariodicose]
    #ATG count excluding start codon
    listadicodoniperfile.extend([aminoacids_sequence.count("M") / boh *100]) # non so cosa sia lol

    listadicodonitotale.append(listadicodoniperfile)

    
dataframe_data = []
temp = []
for index, listaperfile2 in enumerate(listadicodonitotale):
    temp = [FileNames[index], IntergenomicGC[index]]
    for cosa in listaperfile2:
        temp.extend([cosa])
    dataframe_data.append(temp)

df = pd.DataFrame(
    data=dataframe_data,
    columns=[
        "Name",
        "IntergenomicGC",
        "AAA",
        "ATA",
        "ACA",
        "AGA",
        "TTA",
        "TCA",
        "CAA",
        "CTA",
        "CCA",
        "CGA",
        "GAA",
        "GTA",
        "GCA",
        "GGA",
        "AAT",
        "ATT",
        "ACT",
        "AGT",
        "TAT",
        "TTT",
        "TCT",
        "TGT",
        "CAT",
        "CTT",
        "CCT",
        "CGT",
        "GAT",
        "GTT",
        "GCT",
        "GGT",
        "AAC",
        "ATC",
        "ACC",
        "AGC",
        "TAC",
        "TTC",
        "TCC",
        "TGC",
        "CAC",
        "CTC",
        "CCC",
        "CGC",
        "GAC",
        "GTC",
        "GCC",
        "GGC",
        "AAG",
        "ACG",
        "AGG",
        "TTG",
        "TCG",
        "TGG",
        "CAG",
        "CTG",
        "CCG",
        "CGG",
        "GAG",
        "GTG",
        "GCG",
        "GGG",
        "ATG"
    ],
    index=range(1, files_tot+1)
)

folder = "results/"
print(f"Saving data to the {folder} folder")
df.to_csv(f"{folder}data Mio.csv", index=False)