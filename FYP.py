from Bio import SeqIO, Seq
import pandas as pd
from os import walk
import os

# Looking for the files in the /genomes folder
mypath = os.path.abspath(os.getcwd()) + "/genomes"
ext = "embl"
files = []
for (dirpath, dirnames, filenames) in walk(mypath):
    files.extend(filenames)
    break
if not files:  # if files is empty
    raise Exception("Warning there are no filenames, check the path!")
    quit()  # no need to proceed


# Declaring all the variables
FileNames = []
GCcontent = []
IntergenomicGC = []
GC3content = []
PhePer = []
LeuPer = []
IlePer = []
MetPer = []
ValPer = []
SerPer = []
IsoPer = []
ProPer = []
ThrPer = []
AlaPer = []
TyrPer = []
HisPer = []
GlnPer = []
AsnPer = []
LysPer = []
AspPer = []
GluPer = []
CysPer = []
TrpPer = []
ArgPer = []
GlyPer = []

files_tot = len(files)
for index, name in enumerate(files):
    print("\n")
    print(f"Analyzing file {index+1} of {files_tot}: {name}")
    print("\n")
    record = SeqIO.read(f'genomes/{name}', ext)
    FileNames.append(str(name))
    sequence = record.seq

    coding_sequences = []  # It will contain the non pseudogenes (reversed if needed)
    intergenomic_sequence = []  # It will contain the sequences out of the locations bounds as strings
    intergenomic_sequence_list = []  # It will contain the locations bounds
    intergenomic_sequence_string = ""  # Intergenomic sequence as a whole (as a string)

    start = 0

    start = 0 
    for feature in record.features:

        locations = feature.location
        st = int(locations.strand)
        gene_start = int(locations.nofuzzy_start)
        gene_end = int(locations.nofuzzy_end)

        if gene_start > start:
            ig_sequence = sequence[start: gene_start]
            intergenomic_sequence.append(str(ig_sequence))
            intergenomic_sequence_list.append((start, gene_start))

        if gene_end != len(sequence):
            start = gene_end  # Moving the pointer to the end of the gene

        # In bound analysis (non sapevo come chiamarla lol)
        genes = str(sequence[gene_start: gene_end])

        # Reversing complementary strands
        if st == -1:
            sequence_negativestrand = genes[::-1]
            sequence_negativestrand = sequence_negativestrand.replace("A", "b")
            sequence_negativestrand = sequence_negativestrand.replace("T", "A")
            sequence_negativestrand = sequence_negativestrand.replace("b", "T")
            sequence_negativestrand = sequence_negativestrand.replace("C", "d")
            sequence_negativestrand = sequence_negativestrand.replace("G", "C")
            sequence_negativestrand = sequence_negativestrand.replace("d", "G")
            genes = sequence_negativestrand

        # Discarding the pseudogenes
        ok = 0
        if len(genes) % 3 != 0:
            ok = 1
        elif genes[:3] != "ATG" and genes[:3] != "GTG" and genes[:3] != "TTG":
            ok = 2
        elif genes[-3:] != "TAA" and genes[-3:] != "TAG" and genes[-3:] != "TGA":
            ok = 3

        # Saving both the complementary and "normal"(come se chiamano?) genes
        if ok == 0:
            coding_sequences.append(genes)

        # Checking for errors
        check = False
        if check and st == 1 and ok != 0: # Perchè "st == 1"?
            print(locations)
            print(ok)
            print(genes[:3])
            print(genes[-3:])

    # Catching the last part of the sequence after the last gene (if present)
    if gene_end != len(sequence):
        ig_sequence = sequence[gene_end:len(sequence)]
        intergenomic_sequence.append(str(ig_sequence))
        intergenomic_sequence_list.append([gene_end, len(sequence)] )

    intergenomic_sequence_string = "".join(intergenomic_sequence)

    # Genomic GC content
    GC = record.seq.count("G") + record.seq.count("C")
    GC_content_percent = (GC / len(record.seq)) * 100
    print(f"Genomic GC content:  {GC_content_percent:.2f}%")
    GCcontent.append(str(GC_content_percent))

    # Intergenomic GC content
    IgGC = (intergenomic_sequence_string.count("G") +
            intergenomic_sequence_string.count("C"))
    IgGC_content_percent = (IgGC / len(intergenomic_sequence_string)) * 100
    print(f"Intergenomic GC content: {IgGC_content_percent:.2f}%")
    IntergenomicGC.append(str(IgGC_content_percent))

    # GC3 content
    Gcount = 0
    Ccount = 0
    Acount = 0
    Tcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            G3 = codons[2].count("G")
            Gcount += G3
            C3 = codons[2].count("C")
            Ccount += C3
            A3 = codons[2].count("A")
            Acount += A3
            T3 = codons[2].count("T")
            Tcount += T3
    GC3 = ((Gcount + Ccount) / (Gcount + Ccount + Acount + Tcount) * 100)
    print(f'GC3: {(GC3):.2f}%')
    GC3content.append(GC3)  # Perchè avevi fatto str?

    # Amino acids usage
    coding_sequence_whole = "".join(coding_sequences)
    aminoacids_sequence = Seq.Seq(coding_sequence_whole)
    aminoacids_sequence = aminoacids_sequence.translate()  # ah mò lo usi il translate ahhaha
    # number of Valine
    # print(f'Total Valine: {aminoacids_sequence.count("V")}')
    print(f'Percentage of Valine:        {(aminoacids_sequence.count("V") / len(aminoacids_sequence) * 100):.2f}%')
    ValPer.append(aminoacids_sequence.count("V") / len(aminoacids_sequence) * 100)  # Float value
    # number of Phenilalanine
    print(f'Percentage of Phenilalanine: {(aminoacids_sequence.count("F") / len(aminoacids_sequence) * 100):.2f}%')
    PhePer.append(aminoacids_sequence.count("F") / len(aminoacids_sequence) * 100)  # Float value
    # number of Glycine
    print(f'Percentage of Glycine:       {(aminoacids_sequence.count("G") / len(aminoacids_sequence) * 100):.2f}%')
    GlyPer.append(aminoacids_sequence.count("G") / len(aminoacids_sequence) * 100)  # Float value
    # number of Leucine
    print(f'Percentage of Leucine:       {(aminoacids_sequence.count("L") / len(aminoacids_sequence) * 100):.2f}%')
    LeuPer.append(aminoacids_sequence.count("L") / len(aminoacids_sequence) * 100)  # Float value
    # number of Serine
    print(f'Percentage of Serine:        {(aminoacids_sequence.count("S") / len(aminoacids_sequence) * 100):.2f}%')
    SerPer.append(aminoacids_sequence.count("S") / len(aminoacids_sequence) * 100)  # Float value + AGGIUNTA PERCENTUALEEEE!!!
    # number of Isoleucine
    print(f'Percentage of Isoleucine:    {(aminoacids_sequence.count("I") / len(aminoacids_sequence) * 100):.2f}%')
    IsoPer.append(aminoacids_sequence.count("I") / len(aminoacids_sequence) * 100)  # Float value
    # number of Proline
    print(f'Percentage of Proline:       {(aminoacids_sequence.count("P") / len(aminoacids_sequence) * 100):.2f}%')
    ProPer.append(aminoacids_sequence.count("P") / len(aminoacids_sequence) * 100)  # Float value
    # number of Threonine
    print(f'Percentage of Threonine:     {(aminoacids_sequence.count("T") / len(aminoacids_sequence) * 100):.2f}%')
    ThrPer.append(aminoacids_sequence.count("T") / len(aminoacids_sequence) * 100)  # Float value
    # number of Alanine
    print(f'Percentage of Alanine:       {(aminoacids_sequence.count("A") / len(aminoacids_sequence) * 100):.2f}%')
    AlaPer.append(aminoacids_sequence.count("A") / len(aminoacids_sequence) * 100)  # Float value + AGGIUNTA PERCENTUALEEEE!!!
    # number of Tyrosine
    print(f'Percentage of Tyrosine:      {(aminoacids_sequence.count("Y") / len(aminoacids_sequence) * 100):.2f}%')
    TyrPer.append(aminoacids_sequence.count("Y") / len(aminoacids_sequence) * 100)  # Float value
    # number of Histidine
    print(f'Percentage of Histidine:     {(aminoacids_sequence.count("H") / len(aminoacids_sequence) * 100):.2f}%')
    HisPer.append(aminoacids_sequence.count("H") / len(aminoacids_sequence) * 100)  # Float value
    # number of Glutamine
    print(f'Percentage of Glutamine:     {(aminoacids_sequence.count("Q") / len(aminoacids_sequence) * 100):.2f}%')
    GlnPer.append(aminoacids_sequence.count("Q") / len(aminoacids_sequence) * 100)  # Float value
    # number of Asparagine
    print(f'Percentage of Asparagine:    {(aminoacids_sequence.count("N") / len(aminoacids_sequence) * 100):.2f}%')
    AsnPer.append(aminoacids_sequence.count("N") / len(aminoacids_sequence) * 100)  # Float value
    # number of Lysine
    print(f'Percentage of Lysine:        {(aminoacids_sequence.count("K") / len(aminoacids_sequence) * 100):.2f}%')
    LysPer.append(aminoacids_sequence.count("K") / len(aminoacids_sequence) * 100)  # Float value
    # number of Aspartic Acid
    print(f'Percentage of Aspartic Acid: {(aminoacids_sequence.count("D") / len(aminoacids_sequence) * 100):.2f}%')
    AspPer.append(aminoacids_sequence.count("D") / len(aminoacids_sequence) * 100)  # Float value
    # number of Glutamic Acid
    print(f'Percentage of Glutamic Acid: {(aminoacids_sequence.count("E") / len(aminoacids_sequence) * 100):.2f}%')
    GluPer.append(aminoacids_sequence.count("E") / len(aminoacids_sequence) * 100)  # Float value
    # number of Cysteine
    print(f'Percentage of Cysteine:      {(aminoacids_sequence.count("C") / len(aminoacids_sequence) * 100):.2f}%')
    CysPer.append(aminoacids_sequence.count("C") / len(aminoacids_sequence) * 100)  # Float value
    # number of Arginine
    print(f'Percentage of Arginine       {(aminoacids_sequence.count("R") / len(aminoacids_sequence) * 100):.2f}%')
    ArgPer.append(aminoacids_sequence.count("R") / len(aminoacids_sequence) * 100)  # Float value
    # number of Tryptophan
    print(f'Percentage of Arginine       {(aminoacids_sequence.count("W") / len(aminoacids_sequence) * 100):.2f}%')
    TrpPer.append(aminoacids_sequence.count("W") / len(aminoacids_sequence) * 100)  # Float value
    # number of Methionine
    aminoacids_sequence = Seq.Seq(coding_sequence_whole[3:]).translate()
    print(f'Percentage of Methionine:    {(aminoacids_sequence.count("M") / len(aminoacids_sequence) * 100):.2f}%')
    MetPer.append(aminoacids_sequence.count("M") / len(aminoacids_sequence) * 100)  # Float value



    # Codon usage for Arginine + R4 vs R2
    CGAcount = 0
    CGCcount = 0
    CGGcount = 0
    CGTcount = 0
    AGAcount = 0
    AGGcount = 0

    # Codon usage for Valine
    GTAcount = 0
    GTCcount = 0
    GTGcount = 0
    GTTcount = 0

    # Codon usage for Proline
    CCAcount = 0
    CCCcount = 0
    CCGcount = 0
    CCTcount = 0

    # Codon usage for Threonine
    ACAcount = 0
    ACCcount = 0
    ACGcount = 0
    ACTcount = 0

    # Codon usage for Alanine
    GCAcount = 0
    GCCcount = 0
    GCGcount = 0
    GCTcount = 0

    for cds in coding_sequences: # Potrebbero essere uniti al loop di prima volendo 
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]

            # Codon usage for Arginine + R4 vs R2
            CGA = codons.count("CGA")
            CGAcount += CGA
            CGC = codons.count("CGC")
            CGCcount += CGC
            CGG = codons.count("CGG")
            CGGcount += CGG
            CGT = codons.count("CGT")
            CGTcount += CGT
            AGA = codons.count("AGA")
            AGAcount += AGA
            AGG = codons.count("AGG")
            AGGcount += AGG

            # Codon usage for Valine
            GTA = codons.count("GTA")
            GTAcount += GTA
            GTC = codons.count("GTC")
            GTCcount += GTC
            GTG = codons.count("GTG")
            GTGcount += GTG
            GTT = codons.count("GTT")
            GTTcount += GTT

            # Codon usage for Proline
            CCA = codons.count("CCA")
            CCAcount += CCA
            CCC = codons.count("CCC")
            CCCcount += CCC
            CCG = codons.count("CCG")
            CCGcount += CCG
            CCT = codons.count("CCT")
            CCTcount += CCT

            # Codon usage for Threonine
            ACA = codons.count("ACA")
            ACAcount += ACA
            ACC = codons.count("ACC")
            ACCcount += ACC
            ACG = codons.count("ACG")
            ACGcount += ACG
            ACT = codons.count("ACT")
            ACTcount += ACT

            # Codon usage for Alanine
            GCA = codons.count("GCA")
            GCAcount += GCA
            GCC = codons.count("GCC")
            GCCcount += GCC
            GCG = codons.count("GCG")
            GCGcount += GCG
            GCT = codons.count("GCT")
            GCTcount += GCT

    # Codon usage for Arginine + R4 vs R2
    Arginine4Total = CGAcount + CGCcount + CGGcount + CGTcount
    Arginine2Total = AGAcount + CGGcount
    ArginineTotal = Arginine4Total + Arginine2Total
    print(f'Arginine 4 percentage: {(Arginine4Total / (ArginineTotal) * 100):.2f}%')
    print(f'Arginine 2 percentage: {(Arginine2Total / (ArginineTotal) * 100):.2f}%')

    # Codon usage for Valine
    ValineTotal = GTAcount + GTCcount + GTGcount + GTTcount
    print(f'GTA percentage: {(GTAcount / (ValineTotal) * 100):.2f}%')
    print(f'GTC percentage: {(GTCcount / (ValineTotal) * 100):.2f}%')
    print(f'GTG percentage: {(GTGcount / (ValineTotal) * 100):.2f}%')
    print(f'GTT percentage: {(GTTcount / (ValineTotal) * 100):.2f}%')

    # Codon usage for Proline
    ProlineTotal = CCAcount + CCCcount + CCGcount + CCTcount
    print(
        f'Codon usage for Proline:\n CCA percentage: {(CCAcount * 100 / (ProlineTotal)):.2f}% \n CCC percentage: {(CCCcount * 100 / (ProlineTotal)):.2f}% \n CCG percentage: {(CCGcount * 100 / (ProlineTotal)):.2f}% \n CCT percentage: {(CCTcount * 100 / (ProlineTotal)):.2f}%')

    # Codon usage for Threonine
    ThreonineTotal = ACAcount + ACCcount + ACGcount + ACTcount
    print(
        f'Codon usage for Threonine:\n ACA percentage: {(ACAcount * 100 / (ThreonineTotal)):.2f}% \n ACC percentage: {(ACCcount * 100 / (ThreonineTotal)):.2f}% \n ACG percentage: {(ACGcount * 100 / (ThreonineTotal)):.2f}% \n ACT percentage: {(ACTcount * 100 / (ThreonineTotal)):.2f}%')

    # Codon usage for Alanine
    AlanineTotal = GCAcount + GCCcount + GCGcount + GCTcount
    print(
        f'Codon usage for Alanine:\n GCA percentage: {(GCAcount * 100 / (AlanineTotal)):.2f}% \n GCC percentage: {(GCCcount * 100 / (AlanineTotal)):.2f}% \n GCG percentage: {(GCGcount * 100 / (AlanineTotal)):.2f}% \n GCT percentage: {(GCTcount * 100 / (AlanineTotal)):.2f}%')

# Eliminata perchè vuota
# "IlePer": IlePer, 
df = pd.DataFrame(
    {
        "Name": FileNames,
        "GCcontent": GCcontent,
        "IntergenomicGC": IntergenomicGC,
        "GC3content": GC3content,
        "PhePer": PhePer,
        "LeuPer": LeuPer,

        "MetPer": MetPer,
        "ValPer": ValPer,
        "SerPer": SerPer,
        "IsoPer": IsoPer,
        "ProPer": ProPer,
        "ThrPer": ThrPer,
        "AlaPer": AlaPer,
        "TyrPer": TyrPer,
        "HisPer": HisPer,
        "GlnPer": GlnPer,
        "AsnPer": AsnPer,
        "LysPer": LysPer,
        "AspPer": AspPer,
        "GluPer": GluPer,
        "CysPer": CysPer,
        "TrpPer": TrpPer,
        "ArgPer": ArgPer,
        "GlyPer": GlyPer
    }
)

df.head()
df.to_csv("results.csv")
