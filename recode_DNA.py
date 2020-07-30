from Bio.Alphabet import IUPAC
from Bio import Seq, SeqIO, SeqFeature
from CDS import CDS
import copy
import pandas as pd
from reportlab.lib import colors



def set_codon_replacement_list(codon_usage_table):
    xl_file = pd.ExcelFile(codon_usage_table)
    dfs = pd.read_excel(xl_file)
    dfs = dfs[["codon", "reassigned codon"]].dropna()
    codons = dfs["codon"].tolist()
    reasigned_codon = dfs["reassigned codon"].tolist()
    # removes the * from the codons in the file
    codons = ' '.join(codons).replace('*', '').split()
    return dict(zip(codons,reasigned_codon))




# takes in a genbank element, and returns a recoded genbank element, using codon_replacement_list
# , the help is if there are multiple of the same CDS in the gb file
def recode(in_gb="genbank/cp.gb", codon_usage_table =  "excel/codon_reassignment.xlsx", in_help="fasta/genes.fasta"):
    cds = CDS(in_gb, fasta=in_help)
    codon_replacements = set_codon_replacement_list(codon_usage_table)
    genbank = copy.deepcopy(cds.gb)
    # changes all the codons in the CDS_info
    change = []
    new_codon_features = []
    n_codon_changed = 0
    for i in cds.CDS_info:
        codons = i["codons"]
        for codon_index in range(len(codons)):
            str_codon = "".join(codons[codon_index][0])
            if (str_codon in codon_replacements.keys()):
                new_codon = codons[codon_index].copy()
                new_codon[0] = list(codon_replacements[str_codon])
                new_codon_features.append(new_codon_feature(codons[codon_index],new_codon))
                change = change + get_base_change(codons[codon_index], new_codon)
                n_codon_changed += 1
    genbank.features = genbank.features+new_codon_features
    new_seq = list(genbank)
    for i in change:
        new_seq[i[1]] = i[0]

    genbank.seq = Seq.Seq("".join(new_seq), IUPAC.ambiguous_dna)
    print("A total of " + str(n_codon_changed) + " codons were changed!")
    return genbank


# return changes than has to happen to change c1 into c2
def get_base_change(c1, c2):
    change = []
    for i in range(len(c1[0])):
        if (c1[0][i] != c2[0][i]):
            if (c1[2][i]) == "f":
                change.append([c2[0][i], c2[1][i]])
            else:
                change.append([flip_base(c2[0][i]), c2[1][i]])
    return (change)

def flip_base(b):
    b = b.upper()
    if(b == "A"):
        return "T"
    elif(b == "T"):
        return "A"
    elif(b == "G"):
        return "C"
    elif(b == "C"):
        return "G"
    else:
        print(b +" could not be flipped...")
        return b


def new_codon_feature(c1, c2):
    c1_str = "".join(c1[0])
    c2_str = "".join(c2[0])
    type = "modified_base"
    color = colors.burlywood.hexval()[2:]
    label= c1_str + " -> " + c2_str

    if(c1[2][0] == "f"):
        strand = 1
        dir = "RIGHT"
        start = c1[1][0]
        end = c1[1][2]+1
    else:
        strand = -1
        dir = "LEFT"
        start = c1[1][2]
        end = c1[1][0]+1

    # just in case there are condons split between CDS's (end/start of exons)
    if(abs(c1[1][2]-c1[1][0]) != 2):
        start = c1[1][1]
        end = start


    feature = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(start,end), strand=strand,
                                    type=type, qualifiers={"label": [label],
                                                           "note": ['color: #'+color+'; direction: '+dir]})
    feature._set_strand(1)
    return feature

def write_recode_file(recoded_gb, in_gb="genbank/cp.gb", only_cds=True):
    if only_cds:
        new_features = []
        for feature in recoded_gb.features:
            if feature.type == "CDS":
                feature.qualifiers["label"] = ["CR " + feature.qualifiers["label"][0]]
                new_features.append(feature)
            elif feature.type == "modified_base":
                new_features.append(feature)
        recoded_gb.features = new_features
    out_gb_file = in_gb.split("/")[-1]
    out_gb_file = "genbank/cr_"+out_gb_file
    SeqIO.write(recoded_gb, out_gb_file, "genbank")

    return


if __name__ == '__main__':
    infile = "genbank/cp.gb"
    infile_help = "fasta/genes.fasta"
    # codon_usage_table = "excel/codon_usage_tables.xlsx"
    codon_usage_table = "excel/codon_reassignment.xlsx"


    recoded_gb = recode()
    write_recode_file(recoded_gb)
    # cr_cds = CDS("cr_cp.gb",fasta="psbA_cr.fasta")
    # for i in cr_cds.CDS_info:
    #     if i["label"] == "psbA CDS (photosystem II protein D1) 1":
    #         s = ""
    #         for j in i["sequence"]:
    #             s = s+j
    #         print("here comes normal version:")
    #         print(s)
    #     if i["label"] == "psbA CDS (photosystem II protein D1) 2":
    #         s = ""
    #         for j in i["sequence"]:
    #             s = s+j
    #         print("here comes normal version 2:")
    #         print(s)
    # print("hello")
    print("a")