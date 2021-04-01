rom Bio import SeqIO
import pandas
from os import walk
import os

mypath = os.path.abspath(os.getcwd()) + "/genomes"
files = []
for (dirpath, dirnames, filenames) in walk(mypath):
    files.extend(filenames)
    break
ext = "embl"

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


for name in files:
    record = SeqIO.read(f'genomes/{name}', ext)
    FileNames.append(str(name))
    sequence = record.seq
    coding_sequences = []

    start = 0
    intergenomic_sequence = []  # It will contain the sequences out of the locations bounds
    Intergenomic_sequence_list = []  # It will contain the locations bounds
    intergenomic_sequence_string = ""  # intergenomic sequence

    for feature in record.features:
        locations = feature.location
        s = int(locations.nofuzzy_start)
        e = int(locations.nofuzzy_end)
        st = int(locations.strand)

        gene_start = int(locations.nofuzzy_start)

        if gene_start > start:
            ig_sequence = sequence[start:gene_start]
            intergenomic_sequence.append(str(ig_sequence))
            Intergenomic_sequence_list.append((start, gene_start))

        new_start = int(locations.nofuzzy_end)

        if new_start != len(sequence):
            start = new_start

        # Complementary strands
        if st == -1:
            sequence_negativestrand = str(sequence[s:e])
            sequence_negativestrand = sequence_negativestrand[::-1]
            sequence_negativestrand = sequence_negativestrand.replace("A", "b")
            sequence_negativestrand = sequence_negativestrand.replace("T", "A")
            sequence_negativestrand = sequence_negativestrand.replace("b", "T")
            sequence_negativestrand = sequence_negativestrand.replace("C", "d")
            sequence_negativestrand = sequence_negativestrand.replace("G", "C")
            sequence_negativestrand = sequence_negativestrand.replace("d", "G")

            # Discarding the pseudogenes
            ok = 0
            genes = str(sequence_negativestrand)

            if len(genes) % 3 != 0:
                ok = 1
            elif genes[:3] != "ATG" and genes[:3] != "GTG" and genes[:3] != "TTG":
                ok = 2
            elif genes[-3:] != "TAA" and genes[-3:] != "TAG" and genes[-3:] != "TGA":
                ok = 3

            if st == -1 and ok == 0:
                coding_sequences.append(genes)
        ok = 0
        genes = str(sequence[s:e])

        if len(genes) % 3 != 0:
            ok = 1
        elif genes[:3] != "ATG" and genes[:3] != "GTG" and genes[:3] != "TTG":
            ok = 2
        elif genes[-3:] != "TAA" and genes[-3:] != "TAG" and genes[-3:] != "TGA":
            ok = 3
        elif len(genes) == len(sequence):
            ok = 4

        if st == 1 and ok == 0:
            coding_sequences.append(genes)

        # checking for errors
        check = False
        if check and st == 1 and ok != 0:
            print(locations)
            print(ok)
            print(genes[:3])
            print(genes[-3:])

    if start != len(sequence):
        ig_sequence = sequence[start:len(sequence)]
        intergenomic_sequence.append(str(ig_sequence))
        Intergenomic_sequence_list.append((start, gene_start))

    intergenomic_sequence_string = "".join(intergenomic_sequence)

    # Genomic GC content
    GC = record.seq.count("G") + record.seq.count("C")
    GC_content = (GC / len(record.seq)) * 100
    print(f"Genomic GC content:  {GC_content:.2f}%")
    GCcontent.append(str(GC_content))

    # Intergenomic GC content
    IgGC = intergenomic_sequence_string.count("G") + intergenomic_sequence_string.count("C")
    IgGC_content = (IgGC / len(intergenomic_sequence_string)) * 100
    print(f"Intergenomic GC content: {IgGC_content:.2f}%")
    IntergenomicGC.append(str(IgGC_content))

    #GC3 content
    Gcount = 0
    Ccount = 0
    Acount = 0
    Tcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            # while codons[:-3] != "TAA" and codons[:-3] !="TAG" and codons[:-3] !="TGA":
            # Does it get rid of only that codon though?
            G3 = codons[2].count("G")
            Gcount += G3
            C3 = codons[2].count("C")
            Ccount += C3
            A3 = codons[2].count("A")
            Acount += A3
            T3 = codons[2].count("T")
            Tcount += T3
    GC3 = ((Gcount + Ccount) / (Gcount + Ccount + Acount + Tcount) * 100):.2f}%')
    print(f'GC3: {(GC3):.2f}%'
    GC3content.append(str(GC3))

    # Amino acids usage
    coding_sequencewhole = "".join(coding_sequences)
    aminoacids_sequence = Seq(coding_sequencewhole).translate()
    # number of Valine
    # print(f'Total Valine: {aminoacids_sequence.count("V")}')
    print(f'Percentage of Valine:        {(aminoacids_sequence.count("V") / len(aminoacids_sequence) * 100):.2f}%')
    ValPer.append(str(aminoacids_sequence.count("V") / len(aminoacids_sequence) * 100))
    # number of Phenilalanine
    print(f'Percentage of Phenilalanine: {(aminoacids_sequence.count("F") / len(aminoacids_sequence) * 100):.2f}%')
    PhePer.append(str(aminoacids_sequence.count("F") / len(aminoacids_sequence) * 100))
    # number of Glycine
    print(f'Percentage of Glycine:       {(aminoacids_sequence.count("G") / len(aminoacids_sequence) * 100):.2f}%')
    GlyPer.append(str(aminoacids_sequence.count("G") / len(aminoacids_sequence) * 100))
    # number of Leucine
    print(f'Percentage of Leucine:       {(aminoacids_sequence.count("L") / len(aminoacids_sequence) * 100):.2f}%')
    LeuPer.append(str(aminoacids_sequence.count("L") / len(aminoacids_sequence) * 100))
    # number of Serine
    print(f'Percentage of Serine:        {(aminoacids_sequence.count("S") / len(aminoacids_sequence) * 100):.2f}%')
    SerPer.append(str(aminoacids_sequence.count("S") / len(aminoacids_sequence))
    # number of Isoleucine
    print(f'Percentage of Isoleucine:    {(aminoacids_sequence.count("I") / len(aminoacids_sequence) * 100):.2f}%')
    IsoPer.append(str(aminoacids_sequence.count("I") / len(aminoacids_sequence) * 100))
    # number of Proline
    print(f'Percentage of Proline:       {(aminoacids_sequence.count("P") / len(aminoacids_sequence) * 100):.2f}%')
    ProPer.append(str(aminoacids_sequence.count("P") / len(aminoacids_sequence) * 100))
    # number of Threonine
    print(f'Percentage of Threonine:     {(aminoacids_sequence.count("T") / len(aminoacids_sequence) * 100):.2f}%')
    ThrPer.append(str(aminoacids_sequence.count("T") / len(aminoacids_sequence) * 100))
    # number of Alanine
    print(f'Percentage of Alanine:       {(aminoacids_sequence.count("A") / len(aminoacids_sequence) * 100):.2f}%')
    AlaPer.append(str(aminoacids_sequence.count("A") / len(aminoacids_sequence))
    # number of Tyrosine
    print(f'Percentage of Tyrosine:      {(aminoacids_sequence.count("Y") / len(aminoacids_sequence) * 100):.2f}%')
    TyrPer.append(str(aminoacids_sequence.count("Y") / len(aminoacids_sequence) * 100))
    # number of Histidine
    print(f'Percentage of Histidine:     {(aminoacids_sequence.count("H") / len(aminoacids_sequence) * 100):.2f}%')
    HisPer.append(str(aminoacids_sequence.count("H") / len(aminoacids_sequence) * 100))
    # number of Glutamine
    print(f'Percentage of Glutamine:     {(aminoacids_sequence.count("Q") / len(aminoacids_sequence) * 100):.2f}%')
    GlnPer.append(str(aminoacids_sequence.count("Q") / len(aminoacids_sequence) * 100))
    # number of Asparagine
    print(f'Percentage of Asparagine:    {(aminoacids_sequence.count("N") / len(aminoacids_sequence) * 100):.2f}%')
    AsnPer.append(str(aminoacids_sequence.count("N") / len(aminoacids_sequence) * 100))
    # number of Lysine
    print(f'Percentage of Lysine:        {(aminoacids_sequence.count("K") / len(aminoacids_sequence) * 100):.2f}%')
    LysPer.append(str(aminoacids_sequence.count("K") / len(aminoacids_sequence) * 100))
    # number of Aspartic Acid
    print(f'Percentage of Aspartic Acid: {(aminoacids_sequence.count("D") / len(aminoacids_sequence) * 100):.2f}%')
    AspPer.append(str(aminoacids_sequence.count("D") / len(aminoacids_sequence) * 100))
    # number of Glutamic Acid
    print(f'Percentage of Glutamic Acid: {(aminoacids_sequence.count("E") / len(aminoacids_sequence) * 100):.2f}%')
    GluPer.append(str(aminoacids_sequence.count("E") / len(aminoacids_sequence) * 100))
    # number of Cysteine
    print(f'Percentage of Cysteine:      {(aminoacids_sequence.count("C") / len(aminoacids_sequence) * 100):.2f}%')
    CysPer.append(str(aminoacids_sequence.count("C") / len(aminoacids_sequence) * 100))
    # number of Arginine
    print(f'Percentage of Arginine       {(aminoacids_sequence.count("R") / len(aminoacids_sequence) * 100):.2f}%')
    ArgPer.append(str(aminoacids_sequence.count("R") / len(aminoacids_sequence) * 100))
    # number of Tryptophan
    print(f'Percentage of Arginine       {(aminoacids_sequence.count("W") / len(aminoacids_sequence) * 100):.2f}%')
    TrpPer.append(str(aminoacids_sequence.count("W") / len(aminoacids_sequence) * 100))
    # number of Methionine
    aminoacids_sequence = Seq(coding_sequencewhole[3:]).translate()
    print(f'Percentage of Methionine:    {(aminoacids_sequence.count("M") / len(aminoacids_sequence) * 100):.2f}%')
    MetPer.append(str(aminoacids_sequence.count("M") / len(aminoacids_sequence) * 100))

    # Codon usage for Arginine + R4 vs R2
    CGAcount = 0
    CGCcount = 0
    CGGcount = 0
    CGTcount = 0
    AGAcount = 0
    AGGcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            CGA = codons.count("CGA")
            CGAcount += CGA
            CGC = codons.count("CGC")
            CGCcount += CGC
            CGG = codons.count("CGG")
            CGGcount += CGG
            CGT = codons.count("CGT")
            CGTcount += CGT
            Arginine4Total = CGAcount + CGCcount + CGGcount + CGTcount
            AGA = codons.count("AGA")
            AGAcount += AGA
            AGG = codons.count("AGG")
            AGGcount += AGG
            Arginine2Total = AGAcount + CGGcount
            ArginineTotal = Arginine4Total + Arginine2Total
    print(f'Arginine 4 percentage: {(Arginine4Total / (ArginineTotal) * 100):.2f}%')
    print(f'Arginine 2 percentage: {(Arginine2Total / (ArginineTotal) * 100):.2f}%')

    # Codon usage for Valine
    GTAcount = 0
    GTCcount = 0
    GTGcount = 0
    GTTcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            GTA = codons.count("GTA")
            GTAcount += GTA
            GTC = codons.count("GTC")
            GTCcount += GTC
            GTG = codons.count("GTG")
            GTGcount += GTG
            GTT = codons.count("GTT")
            GTTcount += GTT
            ValineTotal = GTAcount + GTCcount + GTGcount + GTTcount
    print(f'GTA percentage: {(GTAcount / (ValineTotal) * 100):.2f}%')
    print(f'GTC percentage: {(GTCcount / (ValineTotal) * 100):.2f}%')
    print(f'GTG percentage: {(GTGcount / (ValineTotal) * 100):.2f}%')
    print(f'GTT percentage: {(GTTcount / (ValineTotal) * 100):.2f}%')

    # Codon usage for Proline
    CCAcount = 0
    CCCcount = 0
    CCGcount = 0
    CCTcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            CCA = codons.count("CCA")
            CCAcount += CCA
            CCC = codons.count("CCC")
            CCCcount += CCC
            CCG = codons.count("CCG")
            CCGcount += CCG
            CCT = codons.count("CCT")
            CCTcount += CCT
            ProlineTotal = CCAcount + CCCcount + CCGcount + CCTcount
    print(
        f'Codon usage for Proline:\n CCA percentage: {(CCAcount * 100 / (ProlineTotal)):.2f}% \n CCC percentage: {(CCCcount * 100 / (ProlineTotal)):.2f}% \n CCG percentage: {(CCGcount * 100 / (ProlineTotal)):.2f}% \n CCT percentage: {(CCTcount * 100 / (ProlineTotal)):.2f}%')

    # Codon usage for Threonine
    ACAcount = 0
    ACCcount = 0
    ACGcount = 0
    ACTcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            ACA = codons.count("ACA")
            ACAcount += ACA
            ACC = codons.count("ACC")
            ACCcount += ACC
            ACG = codons.count("ACG")
            ACGcount += ACG
            ACT = codons.count("ACT")
            ACTcount += ACT
            ThreonineTotal = ACAcount + ACCcount + ACGcount + ACTcount
    print(
        f'Codon usage for Threonine:\n ACA percentage: {(ACAcount * 100 / (ThreonineTotal)):.2f}% \n ACC percentage: {(ACCcount * 100 / (ThreonineTotal)):.2f}% \n ACG percentage: {(ACGcount * 100 / (ThreonineTotal)):.2f}% \n ACT percentage: {(ACTcount * 100 / (ThreonineTotal)):.2f}%')

    # Codon usage for Alanine
    GCAcount = 0
    GCCcount = 0
    GCGcount = 0
    GCTcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            GCA = codons.count("GCA")
            GCAcount += GCA
            GCC = codons.count("GCC")
            GCCcount += GCC
            GCG = codons.count("GCG")
            GCGcount += GCG
            GCT = codons.count("GCT")
            GCTcount += GCT
            AlanineTotal = GCAcount + GCCcount + GCGcount + GCTcount
    print(
        f'Codon usage for Alanine:\n GCA percentage: {(GCAcount * 100 / (AlanineTotal)):.2f}% \n GCC percentage: {(GCCcount * 100 / (AlanineTotal)):.2f}% \n GCG percentage: {(GCGcount * 100 / (AlanineTotal)):.2f}% \n GCT percentage: {(GCTcount * 100 / (AlanineTotal)):.2f}%')
