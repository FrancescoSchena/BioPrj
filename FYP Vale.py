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
IntergenomicGC = []
AAA = []
ATA = []
ACA = []
AGA = []
TTA = []
TCA = []
CAA = []
CTA = []
CCA = []
CGA = []
GAA = []
GTA = []
GCA = []
GGA = []
AAT = []
ATT = []
ACT = []
AGT = []
TAT = []
TTT = []
TCT = []
TGT = []
CAT = []
CTT = []
CCT = []
CGT = []
GAT = []
GTT = []
GCT = []
GGT = []
AAC = []
ATC = []
ACC = []
AGC = [] 
TAC = []
TTC = []
TCC = []
TGC = []
CAC = []
CTC = []
CCC = []
CGC = []
GAC = []
GTC = []
GCC = []
GGC = []
AAG = []
ATG = []
ACG = []
AGG = []
TTG = []
TCG = []
TGG = []
CAG = []
CTG = []
CCG = []
CGG = []
GAG = []
GTG = []
GCG = []
GGG = []


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

    # Intergenomic GC content
    IgGC = intergenomic_sequence_string.count("G") + intergenomic_sequence_string.count("C")
    IgGC_content = (IgGC / len(intergenomic_sequence_string)) * 100
    print(f"Intergenomic GC content: {IgGC_content:.2f}%")
    IntergenomicGC.append(IgGC_content)
    
    coding_sequencewhole = "".join(coding_sequences)
    aminoacids_sequence = Seq(coding_sequencewhole).translate()
    
    
    #Codon counting 
    AAAcount = 0
    AATcount = 0
    AACcount = 0
    AAGcount = 0
    ATAcount = 0
    ATTcount = 0
    ATCcount = 0
    ATGcount = 0
    ACAcount = 0
    ACCcount = 0
    ACGcount = 0
    ACTcount = 0
    AGAcount = 0
    AGTcount = 0
    AGCcount = 0
    AGGcount = 0
    TATcount = 0
    TACcount = 0
    TTAcount = 0
    TTTcount = 0
    TTCcount = 0
    TTGcount = 0
    TCAcount = 0
    TCCcount = 0
    TCGcount = 0
    TCTcount = 0
    TGTcount = 0
    TGCcount = 0
    TGGcount = 0
    CAAcount = 0  
    CATcount = 0
    CACcount = 0
    CAGcount = 0
    CTAcount = 0
    CTCcount = 0
    CTGcount = 0
    CTTcount = 0
    CCAcount = 0
    CCCcount = 0
    CCGcount = 0
    CCTcount = 0
    CGAcount = 0
    CGCcount = 0
    CGGcount = 0
    CGTcount = 0
    GAAcount = 0
    GATcount = 0
    GACcount = 0
    GAGcount = 0
    GTAcount = 0
    GTCcount = 0
    GTGcount = 0
    GTTcount = 0
    GCAcount = 0
    GCCcount = 0
    GCGcount = 0
    GCTcount = 0
    GGTcount = 0
    GGCcount = 0
    GGAcount = 0
    GGGcount = 0
    
    for cds in coding_sequences:
        for i in range(0, len(cds), 3):
            codons = cds[i: i + 3]
            AAA1 = codons.count("AAA")
            AAAcount += AAA1
            AAT1 = codons.count("AAT")
            AATcount += AAT1
            AAC1 = codons.count("AAC")
            AACcount += AAC1
            AAG1 = codons.count("AAG")
            AAGcount += AAG1
            ATA1 = codons.count("ATA")
            ATAcount += ATA1
            ATT1 = codons.count("ATT")
            ATTcount += ATT1
            ATC1 = codons.count("ATC")
            ATCcount += ATC1
            ACA1 = codons.count("ACA")
            ACAcount += ACA1
            ACC1 = codons.count("ACC")
            ACCcount += ACC1
            ACG1 = codons.count("ACG")
            ACGcount += ACG1
            ACT1 = codons.count("ACT")
            ACTcount += ACT1
            AGA1 = codons.count("AGA")
            AGAcount += AGA1
            AGT1 = codons.count("AGT")
            AGTcount += AGT1
            AGC1 = codons.count("AGC")
            AGCcount += AGC1
            AGG1 = codons.count("AGG")
            AGGcount += AGG1
            CAA1 = codons.count("CAA")
            CAAcount += CAA1
            CAT1 = codons.count("CAT")
            CATcount += CAT1
            CAC1 = codons.count("CAC")
            CACcount += CAC1
            CAG1 = codons.count("CAG")
            CAGcount += CAG1
            CCA1 = codons.count("CCA")
            CCAcount += CCA1
            CCC1 = codons.count("CCC")
            CCCcount += CCC1
            CCG1 = codons.count("CCG")
            CCGcount += CCG1
            CCT1 = codons.count("CCT")
            CCTcount += CCT1
            CTA1 = codons.count("CTA")
            CTAcount += CTA1
            CTC1 = codons.count("CTC")
            CTCcount += CTC1
            CTG1 = codons.count("CTG")
            CTGcount += CTG1
            CTT1 = codons.count("CTT")
            CTTcount += CTT1
            CGA1 = codons.count("CGA")
            CGAcount += CGA1
            CGC1 = codons.count("CGC")
            CGCcount += CGC1
            CGG1 = codons.count("CGG")
            CGGcount += CGG1
            CGT1 = codons.count("CGT")
            CGTcount += CGT1
            TAT1 = codons.count("TAT")
            TATcount += TAT1
            TAC1 = codons.count("TAC")
            TACcount += TAC1
            TCA1 = codons.count("TCA")
            TCAcount += TCA1
            TCC1 = codons.count("TCC")
            TCCcount += TCC1
            TCG1 = codons.count("TCG")
            TCGcount += TCG1
            TCT1 = codons.count("TCT")
            TCTcount += TCT1
            TTA1 = codons.count("TTA")
            TTAcount += TTA1
            TTT1 = codons.count("TTT")
            TTTcount += TTT1
            TTC1 = codons.count("TTC")
            TTCcount += TTC1
            TTG1 = codons.count("TTG")
            TTGcount += TTG1
            TGT1 = codons.count("TGT")
            TGTcount += TGT1
            TGC1 = codons.count("TGC")
            TGCcount += TGC1
            TGG1 = codons.count("TGG")
            TGGcount += TGG1
            GAA1 = codons.count("GAA")
            GAAcount += GAA1
            GAT1 = codons.count("GAT")
            GATcount += GAT1
            GAC1 = codons.count("GAC")
            GACcount += GAC1
            GAG1 = codons.count("GAG")
            GAGcount += GAG1
            GTA1 = codons.count("GTA")
            GTAcount += GTA1
            GTC1 = codons.count("GTC")
            GTCcount += GTC1
            GTG1 = codons.count("GTG")
            GTGcount += GTG1
            GTT1 = codons.count("GTT")
            GTTcount += GTT1
            GCA1 = codons.count("GCA")
            GCAcount += GCA1
            GCC1 = codons.count("GCC")
            GCCcount += GCC1
            GCG1 = codons.count("GCG")
            GCGcount += GCG1
            GCT1 = codons.count("GCT")
            GCTcount += GCT1
            GGA1 = codons.count("GGA")
            GGAcount += GGA1
            GGC1 = codons.count("GGC")
            GGCcount += GGC1
            GGG1 = codons.count("GGG")
            GGGcount += GGG1
            GGT1 = codons.count("GGT")
            GGTcount += GGT1
            
    #ATG count excluding start codon
    aminoacids_sequence = Seq(coding_sequencewhole[3:]).translate()
    ATG.append(aminoacids_sequence.count("M") / len(aminoacids_sequence) * 100)
    
    
    AAA.append(AAAcount/ len(aminoacids_sequence)*100)
    AAT.append(AATcount/ len(aminoacids_sequence)*100)
    AAC.append(AACcount/ len(aminoacids_sequence)*100)
    AAG.append(AAGcount/ len(aminoacids_sequence)*100)
    ATA.append(ATAcount/ len(aminoacids_sequence)*100)
    ATT.append(ATTcount/ len(aminoacids_sequence)*100)
    ATC.append(ATCcount/ len(aminoacids_sequence)*100)
    ACA.append(ACAcount/ len(aminoacids_sequence)*100)
    ACC.append(ACCcount/ len(aminoacids_sequence)*100)
    ACG.append(ACGcount/ len(aminoacids_sequence)*100)
    ACT.append(ACTcount/ len(aminoacids_sequence)*100)
    AGA.append(AGAcount/ len(aminoacids_sequence)*100)
    AGT.append(AGTcount/ len(aminoacids_sequence)*100)
    AGC.append(AGCcount/ len(aminoacids_sequence)*100)
    AGG.append(AGGcount/ len(aminoacids_sequence)*100)
    TAT.append(TATcount/ len(aminoacids_sequence)*100)
    TAC.append(TACcount/ len(aminoacids_sequence)*100)
    TTA.append(TTAcount/ len(aminoacids_sequence)*100)
    TTT.append(TTTcount/ len(aminoacids_sequence)*100)
    TTC.append(TTCcount/ len(aminoacids_sequence)*100)
    TTG.append(TTGcount/ len(aminoacids_sequence)*100)
    TCA.append(TCAcount/ len(aminoacids_sequence)*100)
    TCC.append(TCCcount/ len(aminoacids_sequence)*100)
    TCG.append(TCGcount/ len(aminoacids_sequence)*100)
    TCT.append(TCTcount/ len(aminoacids_sequence)*100)
    TGT.append(TGTcount/ len(aminoacids_sequence)*100)
    TGC.append(TGCcount/ len(aminoacids_sequence)*100)
    TGG.append(TGGcount/ len(aminoacids_sequence)*100)
    CAA.append(CAAcount/ len(aminoacids_sequence)*100)
    CAT.append(CATcount/ len(aminoacids_sequence)*100)
    CAC.append(CACcount/ len(aminoacids_sequence)*100)
    CAG.append(CAGcount/ len(aminoacids_sequence)*100)
    CTA.append(CTAcount/ len(aminoacids_sequence)*100)
    CTC.append(CTCcount/ len(aminoacids_sequence)*100)
    CTG.append(CTGcount/ len(aminoacids_sequence)*100)
    CTT.append(CTTcount/ len(aminoacids_sequence)*100)
    CCA.append(CCAcount/ len(aminoacids_sequence)*100)
    CCC.append(CCCcount/ len(aminoacids_sequence)*100)
    CCG.append(CCGcount/ len(aminoacids_sequence)*100)
    CCT.append(CCTcount/ len(aminoacids_sequence)*100)
    CGA.append(CGAcount/ len(aminoacids_sequence)*100)
    CGC.append(CGCcount/ len(aminoacids_sequence)*100)
    CGG.append(CGGcount/ len(aminoacids_sequence)*100)
    CGT.append(CGTcount/ len(aminoacids_sequence)*100)
    GAA.append(GAAcount/ len(aminoacids_sequence)*100)
    GAT.append(GATcount/ len(aminoacids_sequence)*100)
    GAC.append(GACcount/ len(aminoacids_sequence)*100)
    GAG.append(GAGcount/ len(aminoacids_sequence)*100)
    GTA.append(GTAcount/ len(aminoacids_sequence)*100)
    GTC.append(GTCcount/ len(aminoacids_sequence)*100)
    GTG.append(GTGcount/ len(aminoacids_sequence)*100)
    GTT.append(GTTcount/ len(aminoacids_sequence)*100)
    GCA.append(GCAcount/ len(aminoacids_sequence)*100)
    GCC.append(GCCcount/ len(aminoacids_sequence)*100)
    GCG.append(GCGcount/ len(aminoacids_sequence)*100)
    GCT.append(GCTcount/ len(aminoacids_sequence)*100)
    GGT.append(GGTcount/ len(aminoacids_sequence)*100)
    GGA.append(GGAcount/ len(aminoacids_sequence)*100)
    GGG.append(GGGcount/ len(aminoacids_sequence)*100)
    GGC.append(GGCcount/ len(aminoacids_sequence)*100)
    
    
df = pd.DataFrame(
    data=list(zip(
        FileNames,
        IntergenomicGC,
        AAA,
        ATA,
        ACA,
        AGA,
        TTA,
        TCA,
        CAA,
        CTA,
        CCA,
        CGA,
        GAA,
        GTA,
        GCA,
        GGA,
        AAT,
        ATT,
        ACT,
        AGT,
        TAT,
        TTT,
        TCT,
        TGT,
        CAT,
        CTT,
        CCT,
        CGT,
        GAT,
        GTT,
        GCT,
        GGT,
        AAC,
        ATC,
        ACC,
        AGC, 
        TAC,
        TTC,
        TCC,
        TGC,
        CAC,
        CTC,
        CCC,
        CGC,
        GAC,
        GTC,
        GCC,
        GGC,
        AAG,
        ATG,
        ACG,
        AGG,
        TTG,
        TCG,
        TGG,
        CAG,
        CTG,
        CCG,
        CGG,
        GAG,
        GTG,
        GCG,
        GGG)
    ),
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
        "ATG",
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
        "GGG"
    ],
    index=range(1, files_tot+1)
)

folder = "results/"
print(f"Saving data to the {folder} folder")
df.to_csv(f"{folder}data.csv", index=False)