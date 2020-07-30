#################################
#    Kristoffer Holm Hundvin    #
#    Master thesis              #
#################################

from itertools import permutations
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

class CDS:
    def __init__(self, gb, fasta = None):
        """
        initializes the class using a genbank file gb, and if possible a fasta file
        containing some of the genes used. most notabily for the "splicing" of cds'
        are not numbered when there are multiple parts of the gene
        :param gb: genbank file annotated with cds
        :param fasta: fasta file with genes
        """
        if fasta is None:
            self.fasta_file = None
        else:
            self.fasta_file = self.get_fasta_genes(fasta)
        self.gb = SeqIO.read(gb, "genbank")
        self.CDS_info = self.get_CDS()
        self.total_codons = self.get_total_codons()
        self.sequence_list = self.get_sequence_dict()
        self.name = gb

    def get_CDS(self):
        """
        Gets the coding sequences from the genbank file
        :return: a dictionary containg the different genes, and their information like
        index, start position, end postion, codons, sequence, direction and name
        """
        count = 0
        label_count = 1
        gene_count = 1
        CDS_info = []
        names = []
        for feature in self.gb.features:
            # collects the genbank feature of type "CDS"
            if feature.type == "CDS":
                if("label" in feature.qualifiers.keys()):
                    if(len(feature.qualifiers["label"])>1 and type(feature.qualifiers["label"]) == list):
                        feature.qualifiers["label"].sort(key=len,reverse=True)
                    testName = feature.qualifiers["label"]
                    name = feature.qualifiers["label"][0]

                # if using gene instead of label, then splicing will not happen
                elif("gene" in feature.qualifiers.keys()):
                    if (len(feature.qualifiers["gene"]) > 1):
                        feature.qualifiers["gene"].sort(key=len, reverse=True)
                    name = feature.qualifiers["gene"][0]+str(gene_count)
                    feature.qualifiers["label"] = [name]
                    gene_count = gene_count+1
                else:
                    feature.qualifiers["label"] = ["non_labled_"+str(count)]
                    label_count = label_count+1
                    name = feature.qualifiers["label"][0]
                start = feature.location.start.position
                end = feature.location.end.position
                sense = feature.strand

                if (sense == 1):
                    strand = 'f'
                    sequence = self.gb.seq[start:end]
                else:
                    strand = 'r'
                    sequence = self.gb.seq[start:end].reverse_complement()
                    temp = start
                    start = end
                    end = temp

                if (name in names):
                    ex_feature_i = names.index(name)
                    new_info = {'start': start, 'end': end, 'sequence': sequence, 'strand': strand}
                    features = self.combineFeature(CDS_info[ex_feature_i], new_info)
                    CDS_info[ex_feature_i]["start"] = features[0]
                    CDS_info[ex_feature_i]["end"] = features[1]
                    CDS_info[ex_feature_i]["strand"] = features[2]
                    CDS_info[ex_feature_i]["sequence"] = features[3]
                else:
                    # codons return a list of codon sequences, and positions in each codon
                    names.append(name)
                    CDS_info.append({'index': count, 'label': name, 'start': [start], 'end': [end],
                                     'sequence': [sequence], 'strand': [strand]})
                    count += 1

        # orders the sequences in one CDS in the correct order
        for i in CDS_info:
            index = CDS_info.index(i)
            if (len(i["sequence"]) > 1):
                seq_with_pos = self.reorderSeq(i["sequence"], i["start"], i["end"], i["strand"])
                start = []
                end = []
                seq = []
                for j in seq_with_pos:
                    seq.append(j[0])
                    start.append(j[1])
                    end.append(j[2])
                CDS_info[index]["start"] = start
                CDS_info[index]["end"] = end
                CDS_info[index]["sequence"] = seq
            CDS_info[index]["codons"] = self.get_codons(CDS_info[index]["sequence"], CDS_info[index]["start"],
                                                   CDS_info[index]["strand"])

        # sorts in the order of appearing CDS in the genbank file, and not alphabetically of the features
        CDS_info = sorted(CDS_info, key=lambda k: k['start'][0])
        return CDS_info

    def get_fasta_genes(self, fasta_file):
        genes = []
        for i in SeqIO.parse(fasta_file, "fasta"):
            genes.append(i.upper())
        return(genes)

    def get_total_codons(self):
        tot_cod = 0
        for i in self.CDS_info:
            tot_cod += len(i["codons"])
        return tot_cod

    def combineFeature(self, f1, f2):
        start = f1["start"]
        start.append(f2["start"])
        end = f1["end"]
        end.append(f2["end"])
        strand = f1["strand"]
        strand.append(f2["strand"])
        sequence = f1["sequence"]
        sequence.append(f2["sequence"])
        return [start, end, strand, sequence]

    # seqList is on the format of [list of possible sequences, list of start for those sequenecs, list of end]
    def reorderSeq(self, seqList, startList, endList, strandList):
        fasta_seqs = []
        if(self.fasta_file != None):
            for i in self.fasta_file:
                fasta_seqs.append(i.seq.upper())

        seq_start_end_list = []
        for i in range(len(seqList)):
            seq_start_end_list.append([seqList[i], startList[i], endList[i], strandList[i]])
        perms = permutations(seq_start_end_list)

        possible_sequences = []
        for i in perms:
            start = Bio.Seq.translate(i[0][0][0:3])
            stop = Bio.Seq.translate(i[-1][0][-3:])
            if (start == "M" and stop == "*"):
                possible_sequence = "".join([str(x[0]) for x in i])
                aa_seq = Bio.Seq.translate(possible_sequence)
                if (aa_seq.count("*") == 1):
                    possible_sequences.append(i)
                    if (possible_sequence in fasta_seqs):
                        return i

        if (len(possible_sequences) == 1):
            return possible_sequences[0]

        elif (len(possible_sequences) == 0):
            # print((seqList[0]))
            print("no possible sequences! ooups...")
            print("there are possibly more than one of the same CDS, with the same name")
            return seq_start_end_list


        elif (len(possible_sequences) > 1):
            print("there are more than one possible combinations of sequences, consider adding a fasta file with the "
                  "correct sequence")
            print("what the hell")

            if(self.fasta_file != None):
                list_of_fasta_AA =[]
                for i in fasta_seqs:
                    list_of_fasta_AA.append(Bio.Seq.translate(str(i)))

                for i in possible_sequences:
                    current_seq = ""
                    for j in i:
                        current_seq = current_seq+j[0]
                        current_AA_seq = Bio.Seq.translate(current_seq)
                    if(current_AA_seq in list_of_fasta_AA):
                        return i

            return None

    # get the codons of the combindes sequences in the seq list, using the start list and the end list
    def get_codons(self, seq, start, strand):
        if (type(seq) != list):
            seq = [seq]
            start = [start]
            strand = [strand]
        codons = []
        triplet = []
        position = []
        base_strand = []
        for i in range(len(seq)):
            for base_index in range(len(seq[i])):

                if (len(triplet) < 3):
                    triplet.append(seq[i][base_index])
                    if (strand[i] == "f"):
                        position.append(start[i] + base_index)
                        base_strand.append("f")
                    else:
                        position.append(start[i] - base_index - 1)
                        base_strand.append("r")

                if (len(triplet) == 3):
                    codons.append([triplet, position, base_strand])
                    triplet = []
                    position = []
                    base_strand = []
        return codons


    def get_sequence_dict(self):
        sl = {}
        for i in self.CDS_info:
            sl[i["label"]] = i["sequence"]
        return sl

    def get_codon_freq(self, codon, cd_info = None):
        """
        :param codon: a three letter string that will be measured freq of
        :param cd_info: an optional parameter that, must be a coding sequence with codons, if included the frequency
        will be measured only in this CDS
        :return: return the frequency of codon in either entire genome, or in the selected CDS per 1000 codons
        """
        counter = 0
        if(cd_info == None):
            for i in self.CDS_info:
                for j in i["codons"]:
                    if ("".join(j[0])) == codon:
                        counter += 1
            div = 1000/self.total_codons
            freq = counter*div
            #freq = counter / self.total_codons
        else:
            for i in cd_info["codons"]:
                if ("".join(i[0])) == codon:
                    counter += 1
            div = 1000 / len(cd_info["codons"])
            freq = counter*div
            #freq = counter / len(cd_info["codons"])
        return freq, counter/self.total_codons


    def get_codon_pair_freq(self, codon1, codon2, cd_info = None):
        counter = 0
        if(cd_info == None):
            for i in self.CDS_info:
                prev = False
                for j in i["codons"]:
                    current = "".join(j[0])
                    if prev:
                        if current == codon2:
                            counter += 1
                    if current == codon1:
                        prev = True
                    else:
                        prev = False
            freq = counter / self.total_codons
        else:
            prev = False
            for i in cd_info["codons"]:
                current = "".join(i[0])
                if prev:
                    if current == codon2:
                        counter += 1
                if current == codon1:
                    prev = True
                else:
                    prev = False
            freq = counter / len(cd_info["codons"])
        return freq

    def get_start_stop_for_extract(self, startGene, endGene):
        go = False
        for i in self.CDS_info:
            if(i["label"][0:3] == "CR " and startGene[0:3] != "CR "):
                startGene = "CR "+startGene
                endGene = "CR "+endGene

            if (i["label"][0:len(startGene)] == startGene and not go):
                start = min(i["start"][0],i["end"][-1])
                go = True
            if(i["label"][0:len(endGene)] == endGene and go):
                end = max(i["end"][-1],i["start"][0])
                return start,end
        return None,None

    def make_genbank_extract(self, startGene,endGene = None):
        if endGene == None:
            endGene = startGene
            filename = "genbank/extractions/ex_" + startGene.replace(" ","")+".gb"
        else:
            filename = "genbank/extractions/ex_" + startGene.replace(" ", "") + "_to_" + endGene.replace(" ", "") + ".gb"
        start,stop = self.get_start_stop_for_extract(startGene, endGene)
        if(start == None):
            print("fuck")
            return

        SeqIO.write(self.gb[start:stop], filename, "genbank")

    # creates a sheet with the codon frequencies
    def get_codon_frequency(self):
        codons = {}
        for i in cds_info:
            for j in i["codons"]:
                cod = "".join(j[0])
                if(cod not in codons.keys()):
                    codons[cod], hold = self.get_codon_freq(cod)

        df = pd.DataFrame(data=codons, index=[0])
        df = (df.T)
        print(df)
        df.to_excel('excel/codon_frequencies.xlsx')
        return codons



# This is for testing
if __name__ == '__main__':

    #cds = CDS("genbank\eColiK-12MG1655.gb")#, "fasta\genes.fasta")
    cds = CDS("genbank\cp.gb", "fasta\genes.fasta")
    #print(cds.get_codon_pair_freq("CTT","CTT", cds.CDS_info[0]))
    #gb = SeqIO.read("test-sequence.gb", "genbank")
    #feature = gb.features[1]
    #cd = CD(feature,gb)
    #cds.draw_full_gene_map()
    # cds.draw_genes_w_DFV()
    cds_info = cds.get_CDS()
    shortest_seq = None
    longest_seq = None
    print(len(cds_info))
    for i in cds_info:
        seq = ""
        for j in i["sequence"]:
            seq = seq+j
        if (shortest_seq == None):
            shortest_seq = seq
        if (longest_seq == None):
            longest_seq = seq
        if(len(seq)<(len(shortest_seq))):
            shortest_seq = seq
        if (len(seq) > (len(longest_seq))):
            longest_seq = seq
       # print(len(seq))
        # print(i["label"])
    # print(len(shortest_seq))
    # print(len(longest_seq))
    for i in cds_info:
        if("rbcL" in i["label"]):
            print(cds.get_codon_freq("GGT", i))

    print(cds.get_codon_frequency())



    # a check to make the figures again if i want
    if True:
        cds.make_genbank_extract("psbH", "psbB")

        cds.make_genbank_extract("psbB")
        cds.make_genbank_extract("psbT")
        cds.make_genbank_extract("psbN")
        cds.make_genbank_extract("psbH")

        cds.make_genbank_extract("rbcL", "psbI")

        cds.make_genbank_extract("rbcL")
        cds.make_genbank_extract("atpA")
        cds.make_genbank_extract("psbI")
