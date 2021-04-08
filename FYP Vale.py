from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

from os import walk
import os
mypath = os.path.abspath(os.getcwd()) + "/filtered_genomes"

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
Leu4Per = []
Leu2Per = []
IlePer = []
MetPer = []
ValPer = []
Ser4Per = []
Ser2Per = []
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
Arg2Per = []
Arg4Per = []
GlyPer = []
PheGC = []
ValGC = []
Ser4GC =[]
Ser2GC =[]
ProGC = []
ThrGC = []
AlaGC = []
TyrGC = []
HisGC = []
GlnGC = []
AsnGC = []
LysGC = []
AspGC = []
GluGC = []
CysGC = []
Arg4GC = []
Arg2GC = []
Leu4GC = []
Leu2GC = []
GlyGC = []

files_tot = len(files)
for index, name in enumerate(files):
    print("\n")
    print(f"Analyzing file {index+1} of {files_tot}: {name}")
    print("\n")
    record = SeqIO.read(f'filtered_genomes/{name}', ext)
    sequence = record.seq
    FileNames.append(str(name))
    
    start = 0
    coding_sequences = []
    intergenomic_sequence = []  # It will contain the sequences out of the locations bounds
    intergenomic_sequence_list = []  # It will contain the locations bounds
    intergenomic_sequence_string = ""  # intergenomic sequence
    
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

        #Reversing complementary strands
        if st == -1:
            sequence_negativestrand = str(sequence[s:e])
            sequence_negativestrand= sequence_negativestrand[::-1]
            sequence_negativestrand= sequence_negativestrand.replace("A","b")
            sequence_negativestrand= sequence_negativestrand.replace("T","A")
            sequence_negativestrand= sequence_negativestrand.replace("b","T")
            sequence_negativestrand= sequence_negativestrand.replace("C","d")
            sequence_negativestrand= sequence_negativestrand.replace("G","C")
            sequence_negativestrand= sequence_negativestrand.replace("d","G")

        #Discarding the pseudogenes

            ok = 0
            genes = str(sequence_negativestrand)


            if len(genes)%3 != 0:
                ok= 1
            elif genes[:3] != "ATG" and genes[:3] != "GTG" and genes[:3] != "TTG":
                ok= 2
            elif genes[-3:] != "TAA" and genes[-3:] != "TAG" and genes[-3:] != "TGA":
                ok= 3


            if st == -1 and ok == 0:
                coding_sequences.append(genes)
        
        ok = 0
        genes = str(sequence[s:e])
    
        if len(genes)%3 != 0:
            ok= 1
        elif genes[:3] != "ATG" and genes[:3] != "GTG" and genes[:3] != "TTG":
            ok= 2
        elif genes[-3:] != "TAA" and genes[-3:] != "TAG" and genes[-3:] != "TGA":
            ok= 3
        elif len(genes) == len(sequence):
            ok= 4

        if st == 1 and ok == 0:
            coding_sequences.append(genes)

        #checking for errors
        check = False
        if check and st == 1 and ok != 0:
            print(locations)
            print(ok)
            print(genes[:3])
            print(genes[-3:])
            
    if start != len(sequence):
        ig_sequence = sequence[start:len(sequence)]
        intergenomic_sequence.append(str(ig_sequence))
        intergenomic_sequence_list.append((start, gene_start))

    intergenomic_sequence_string = "".join(intergenomic_sequence)

    # Genomic GC content
    GC = record.seq.count("G") + record.seq.count("C")
    GC_content = (GC / len(record.seq)) * 100
    print(f"Genomic GC content:  {GC_content:.2f}%")
    GCcontent.append(GC_content)

    # Intergenomic GC content
    IgGC = intergenomic_sequence_string.count("G") + intergenomic_sequence_string.count("C")
    IgGC_content = (IgGC / len(intergenomic_sequence_string)) * 100
    print(f"Intergenomic GC content: {IgGC_content:.2f}%")
    IntergenomicGC.append(IgGC_content)
    
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
    GC3 = (Gcount + Ccount) / (Gcount + Ccount + Acount + Tcount) * 100
    print(f'GC3: {(GC3):.2f}%')
    GC3content.append(GC3)
          
    # Amino acids usage
    coding_sequencewhole = "".join(coding_sequences)
    aminoacids_sequence = Seq(coding_sequencewhole).translate()
    # number of Valine
    print(f'Percentage of Valine:        {(aminoacids_sequence.count("V") / len(aminoacids_sequence) * 100):.2f}%')
    ValPer.append(aminoacids_sequence.count("V") / len(aminoacids_sequence) * 100)
    # number of Phenilalanine
    print(f'Percentage of Phenilalanine: {(aminoacids_sequence.count("F") / len(aminoacids_sequence) * 100):.2f}%')
    PhePer.append(aminoacids_sequence.count("F") / len(aminoacids_sequence) * 100)
    # number of Glycine
    print(f'Percentage of Glycine:       {(aminoacids_sequence.count("G") / len(aminoacids_sequence) * 100):.2f}%')
    GlyPer.append(aminoacids_sequence.count("G") / len(aminoacids_sequence) * 100)
    # number of Isoleucine
    print(f'Percentage of Isoleucine:    {(aminoacids_sequence.count("I") / len(aminoacids_sequence) * 100):.2f}%')
    IlePer.append(aminoacids_sequence.count("I") / len(aminoacids_sequence) * 100)
    # number of Proline
    print(f'Percentage of Proline:       {(aminoacids_sequence.count("P") / len(aminoacids_sequence) * 100):.2f}%')
    ProPer.append(aminoacids_sequence.count("P") / len(aminoacids_sequence) * 100)
    # number of Threonine
    print(f'Percentage of Threonine:     {(aminoacids_sequence.count("T") / len(aminoacids_sequence) * 100):.2f}%')
    ThrPer.append(aminoacids_sequence.count("T") / len(aminoacids_sequence) * 100)
    # number of Alanine
    print(f'Percentage of Alanine:       {(aminoacids_sequence.count("A") / len(aminoacids_sequence) * 100):.2f}%')
    AlaPer.append(aminoacids_sequence.count("A") / len(aminoacids_sequence) * 100)
    # number of Tyrosine
    print(f'Percentage of Tyrosine:      {(aminoacids_sequence.count("Y") / len(aminoacids_sequence) * 100):.2f}%')
    TyrPer.append(aminoacids_sequence.count("Y") / len(aminoacids_sequence) * 100)
    # number of Histidine
    print(f'Percentage of Histidine:     {(aminoacids_sequence.count("H") / len(aminoacids_sequence) * 100):.2f}%')
    HisPer.append(aminoacids_sequence.count("H") / len(aminoacids_sequence) * 100)
    # number of Glutamine
    print(f'Percentage of Glutamine:     {(aminoacids_sequence.count("Q") / len(aminoacids_sequence) * 100):.2f}%')
    GlnPer.append(aminoacids_sequence.count("Q") / len(aminoacids_sequence) * 100)
    # number of Asparagine
    print(f'Percentage of Asparagine:    {(aminoacids_sequence.count("N") / len(aminoacids_sequence) * 100):.2f}%')
    AsnPer.append(aminoacids_sequence.count("N") / len(aminoacids_sequence) * 100)
    # number of Lysine
    print(f'Percentage of Lysine:        {(aminoacids_sequence.count("K") / len(aminoacids_sequence) * 100):.2f}%')
    LysPer.append(aminoacids_sequence.count("K") / len(aminoacids_sequence) * 100)
    # number of Aspartic Acid
    print(f'Percentage of Aspartic Acid: {(aminoacids_sequence.count("D") / len(aminoacids_sequence) * 100):.2f}%')
    AspPer.append(aminoacids_sequence.count("D") / len(aminoacids_sequence) * 100)
    # number of Glutamic Acid
    print(f'Percentage of Glutamic Acid: {(aminoacids_sequence.count("E") / len(aminoacids_sequence) * 100):.2f}%')
    GluPer.append(aminoacids_sequence.count("E") / len(aminoacids_sequence) * 100)
    # number of Cysteine
    print(f'Percentage of Cysteine:      {(aminoacids_sequence.count("C") / len(aminoacids_sequence) * 100):.2f}%')
    CysPer.append(aminoacids_sequence.count("C") / len(aminoacids_sequence) * 100)
    # number of Tryptophan
    print(f'Percentage of Tryptophan     {(aminoacids_sequence.count("W") / len(aminoacids_sequence) * 100):.2f}%')
    TrpPer.append(aminoacids_sequence.count("W") / len(aminoacids_sequence) * 100)
    # number of Methionine
    aminoacids_sequence = Seq(coding_sequencewhole[3:]).translate()
    print(f'Percentage of Methionine:    {(aminoacids_sequence.count("M") / len(aminoacids_sequence) * 100):.2f}%')
    MetPer.append(aminoacids_sequence.count("M") / len(aminoacids_sequence) * 100)
    

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
    ValGC.append((GTGcount + GTCcount) / (ValineTotal) * 100)

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
    ProGC.append((CCCcount + CCGcount) * 100 / (ProlineTotal))

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
    ThrGC.append((ACCcount + ACGcount)* 100 / (ThreonineTotal))

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
    AlaGC.append((GCGcount + GCC) * 100 / (AlanineTotal))
    
    # Codon usage for Phenylalanine
    TTTcount = 0
    TTCcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            TTT = codons.count("TTT")
            TTTcount += TTT
            TTC = codons.count("TTC")
            TTCcount += TTC
    PhenylalanineTotal = TTTcount + TTCcount
    PheGC.append(TTCcount * 100 / PhenylalanineTotal)
    
    # Codon usage for Tyrosine
    TATcount = 0
    TACcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            TAT = codons.count("TAT")
            TATcount += TAT
            TAC = codons.count("TAC")
            TACcount += TAC
    TyrosineTotal = TATcount + TACcount
    TyrGC.append(TACcount * 100 / TyrosineTotal)
    
    # Codon usage for Histidine
    CATcount = 0
    CACcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            CAT = codons.count("CAT")
            CATcount += CAT
            CAC = codons.count("CAC")
            CACcount += CAC
            HistidineTotal = CATcount + CACcount
    HisGC.append(CACcount * 100 / HistidineTotal)
    
    # Codon usage for Glutamine
    CAAcount = 0
    CAGcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            CAA = codons.count("CAA")
            CAAcount += CAA
            CAG = codons.count("CAG")
            CAGcount += CAG
    GlutamineTotal = CAAcount + CAGcount
    GlnGC.append(CAGcount * 100 / GlutamineTotal)    
    
    # Codon usage for Asparagine
    AATcount = 0
    AACcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            AAT = codons.count("AAT")
            AATcount += AAT
            AAC = codons.count("AAC")
            AACcount += AAC
    AsparagineTotal = AATcount + AACcount
    AsnGC.append(AACcount * 100 / AsparagineTotal)
    
    # Codon usage for Lysine
    AAAcount = 0
    AAGcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            AAA = codons.count("AAA")
            AAAcount += AAA
            AAG = codons.count("AAG")
            AAGcount += AAG
    LysineTotal = AAAcount + AAGcount
    LysGC.append(AAGcount * 100 / LysineTotal)
    
    # Codon usage for Aspartic Acid
    GATcount = 0
    GACcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            GAT = codons.count("GAT")
            GATcount += GAT
            GAC = codons.count("GAC")
            GACcount += GAC
    AsparticTotal = GATcount + GACcount
    AspGC.append(GACcount * 100 / AsparticTotal)
    
    # Codon usage for Glutamic Acid
    GAAcount = 0
    GAGcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            GAA = codons.count("GAA")
            GAAcount += GAA
            GAG = codons.count("GAG")
            GAGcount += GAG
    GlutamicTotal = GAAcount + GAGcount
    GluGC.append(GAGcount * 100 / GlutamicTotal)
    
    # Codon usage for Cysteine
    TGTcount = 0
    TGCcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            TGT = codons.count("TGT")
            TGTcount += TGT
            TGC = codons.count("TGC")
            TGCcount += TGC
    CysteineTotal = TGTcount + TGCcount
    CysGC.append(TGCcount * 100 / CysteineTotal)
    
    # Codon usage for Glycine
    GGTcount = 0
    GGCcount = 0
    GGAcount = 0
    GGGcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            GGA = codons.count("GGA")
            GGAcount += GGA
            GGC = codons.count("GGC")
            GGCcount += GGC
            GGG = codons.count("GGG")
            GGGcount += GGG
            GGT = codons.count("GGT")
            GGTcount += GGT
    GlycineTotal = GGAcount + GGCcount + GGGcount + GGTcount
    GlyGC.append((GGGcount + GGC) * 100 / (GlycineTotal))
    
    
    # Codon usage for Arginine, Leucine and Serine
    CGAcount = 0
    CGCcount = 0
    CGGcount = 0
    CGTcount = 0
    AGAcount = 0
    AGGcount = 0
    CTAcount = 0
    CTCcount = 0
    CTGcount = 0
    CTTcount = 0
    TTAcount = 0
    TTGcount = 0
    TCAcount = 0
    TCCcount = 0
    TCGcount = 0
    TCTcount = 0
    AGTcount = 0
    AGCcount = 0
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            #Codon usage for Arginine
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
            #Codon usage for Leucine
            CTA = codons.count("CTA")
            CTAcount += CTA
            CTC = codons.count("CTC")
            CTCcount += CTC
            CTG = codons.count("CTG")
            CTGcount += CTG
            CTT = codons.count("CTT")
            CTTcount += CTT
            TTA = codons.count("TTA")
            TTAcount += TTA
            TTG = codons.count("TTG")
            TTGcount += TTG
            #Codon usage for Serine
            TCA = codons.count("TCA")
            TCAcount += TCA
            TCC = codons.count("TCC")
            TCCcount += TCC
            TCG = codons.count("TCG")
            TCGcount += TCG
            TCT = codons.count("TCT")
            TCTcount += TCT
            AGT = codons.count("AGT")
            AGTcount += AGT
            AGC = codons.count("AGC")
            AGCcount += AGC
    Arginine4Total = CGGcount + CGAcount + CGCcount + CGTcount
    Arginine2Total = AGAcount + AGGcount
    Arg4Per.append(Arginine4Total/len(aminoacids_sequence)*100)
    Arg2Per.append(Arginine2Total/len(aminoacids_sequence)*100)
    Leucine4Total = CTGcount + CTAcount + CTCcount + CTTcount
    Leucine2Total = TTAcount + TTGcount
    Leu4Per.append(Leucine4Total/len(aminoacids_sequence)*100)
    Leu2Per.append(Leucine2Total/len(aminoacids_sequence)*100)
    Serine4Total = TCGcount + TCAcount + TCCcount + TCTcount
    Serine2Total = AGTcount + AGCcount
    Ser4Per.append(Serine4Total/len(aminoacids_sequence)*100)
    Ser2Per.append(Serine2Total/len(aminoacids_sequence)*100)
    
    
    
df = pd.DataFrame(
    data=list(zip(
        FileNames,
        GCcontent,
        IntergenomicGC,
        GC3content,
        PhePer,
        Leu2Per,
        Leu4Per,
        IlePer,
        MetPer,
        ValPer,
        Ser2Per,
        Ser4Per,
        ProPer,
        ThrPer,
        AlaPer,
        TyrPer,
        HisPer,
        GlnPer,
        AsnPer,
        LysPer,
        AspPer,
        GluPer,
        CysPer,
        TrpPer,
        Arg2Per,
        Arg4Per,
        GlyPer,
        PheGC,
        Leu2GC,
        Leu4GC,
        ValGC,
        Ser2GC,
        Ser4GC,
        ProGC,
        ThrGC,
        AlaGC,
        TyrGC,
        HisGC,
        GlnGC,
        AsnGC,
        LysGC,
        AspGC,
        GluGC,
        CysGC,
        Arg2GC,
        Arg4GC,
        GlyGC)
    ),
    columns=[
        "Name",
        "GCcontent",
        "IntergenomicGC",
        "GC3content",
        "PhePer",
        "Leu2Per",
        "Leu4Per",
        "IlePer",
        "MetPer",
        "ValPer",
        "Ser2Per",
        "Ser4Per",
        "ProPer",
        "ThrPer",
        "AlaPer",
        "TyrPer",
        "HisPer",
        "GlnPer",
        "AsnPer",
        "LysPer",
        "AspPer",
        "GluPer",
        "CysPer",
        "TrpPer",
        "Arg2Per",
        "Arg4Per",
        "GlyPer",
        "PheGC",
        "Leu2GC",
        "Leu4GC",
        "ValGC",
        "Ser2GC",
        "Ser4GC",
        "ProGC",
        "ThrGC",
        "AlaGC",
        "TyrGC",
        "HisGC",
        "GlnGC",
        "AsnGC",
        "LysGC",
        "AspGC",
        "GluGC",
        "CysGC",
        "Arg2GC",
        "Arg4GC",
        "GlyGC"
    ],
    index=range(1, files_tot+1)
)

folder = "results/"
print(f"Saving data to the {folder} folder")
df.to_csv(f"{folder}data.csv", index=False)