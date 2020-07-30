from CDS import CDS
from Bio import SeqIO
from Bio.SeqUtils import GC
import pandas as pd
import numpy as np


def make_NN_sequence(cds,gb, before, after):
    start = cds["start"][0]
    strand = cds["strand"][0]
    label = cds["label"]

    if strand == "f":
        new_seq_start = (start-before)%len(gb)
        new_seq_end = (start+after+3)%len(gb)
        if(start < before):
            new_seq = gb[new_seq_start:]+gb[:new_seq_end]
        else:
            new_seq = gb[new_seq_start:new_seq_end]
    else:
        new_seq_start = (start+before)%len(gb)
        new_seq_end = (start-(after+3))%len(gb)
        if(start > len(gb)-(before+4)):
            new_seq = gb[new_seq_end:len(gb)]+gb[0:new_seq_start].reverse_complement()
        else:
            new_seq = gb[new_seq_end:new_seq_start].reverse_complement()
    id_end = label.find("(")
    if(id_end != -1):
        id = label[:id_end]
        desc = label[id_end:]
    else:
        id = label
        desc = ""
    new_seq.id = id + desc
    gc = get_NN_GC_content(new_seq, before, after)
    new_seq.description = desc + " " +str(gc)
    new_seq.gc = gc
    return new_seq

#returns the GC content of the nucs N before ATG, the nucs N after and the total
def get_NN_GC_content(seqNN, before, after):
    start = seqNN.seq[:before+1]
    end = seqNN.seq[before+3:]
    gc = {}
    gc["GC before ATG"] = round(GC(start),3)
    gc["GC after ATG"]= round(GC(end),3)
    gc["total GC"] = round(GC(seqNN.seq),3)
    return gc

def write_NN_fasta(cds, before, after = 0, xlsx = False):
    new_sequences = []
    CDS_info = cds.CDS_info
    outname = cds.name
    outname = outname.rsplit("/")[-1]
    outname =outname.rsplit(".")[0]
    outfile = "_NN_seqs_" + outname + ".fasta"
    if(xlsx):
        outfile_xlsx =  "_NN_seqs_" + outname + ".xlsx"
        total_gc = 0
        ids = []
        gcs = []
    for i in CDS_info:
        new_sequence = make_NN_sequence(i, cds.gb, before, after)
        new_sequences.append(new_sequence)
        if (xlsx):
            current_gc = new_sequence.gc["GC before ATG"]
            total_gc = total_gc + current_gc
            ids.append(new_sequence.id)
            gcs.append(current_gc)
        if(after == 0):
            filename = "fasta/"+ str(before) + outfile
            if (xlsx):
                filename_xlsx = "excel/"+ str(before) + outfile_xlsx
        else:
            filename = "fasta/"+ str(before) + "_" + str(after) + outfile
            if (xlsx):
                filename_xlsx = "excel/"+ str(before) + "_" + str(after) + outfile_xlsx
        SeqIO.write(new_sequences, filename, "fasta")
    if(xlsx):
        columns = ["IDs", "%GC","Varriance", "Outlier"]
        avg_gc = total_gc / len(CDS_info)
        varriance = []
        outlier = []
        q75, q25 = np.percentile(gcs, [75, 25])
        iqr = q75 - q25
        for i in gcs:
            varriance.append((i-avg_gc)**2)
            if(i<= avg_gc+(iqr*1.5) and i>= avg_gc-(iqr*1.5)):
                outlier.append("no")
            else:
                outlier.append("yes")
        print(q25, q75)
        rows = list(zip(ids, gcs,varriance,outlier))
        #rows = [ids,gcs]

        data = pd.DataFrame(rows, columns=columns)
        data.to_excel(filename_xlsx)


if __name__ == '__main__':
    infile = "genbank/cp.gb"
    infile_helper = "fasta/genes.fasta"

    #list of all the CDSs in the genbank file
    cds = CDS(infile, infile_helper)

    # write_NN_fasta(cds, 25)
    # write_NN_fasta(cds, 50, 50)
    write_NN_fasta(cds,200, xlsx=True)





