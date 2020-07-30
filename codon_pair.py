from CDS import CDS
import os.path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def set_codon_lists():
    bases = "TCAG"
    return [a + b + c for a in bases for b in bases for c in bases]


#def make_excel_file(codon_list, in_gb_file, cds):
def make_excel_file( in_gb_file, cds):
    # sets up a lsit with all 64 codons
    codon_list = set_codon_lists()
    in_gb_nofolder = in_gb_file.split("\\")[-1]
    excel_out_name ="excel\\codon_pair_" + in_gb_nofolder.split(".")[0] + ".xlsx"
    pair_values_name = "pair values"
    expected_values_name = "expected values"
    diff_name = "diff"
    return_file = {}
    sheet_names = {}
    if os.path.isfile(excel_out_name):
        print("found " + excel_out_name)
        return_file = pd.read_excel(excel_out_name, sheet_name=None, index_col=0)
    else:
        print("could not find file")
        print("creating " + excel_out_name)
        rows = codon_list
        columns = codon_list
        dict_vals_pairs = {}
        dict_vals_expected = {}
        for row in rows:
            row_list_pairs = []
            row_list_expected = []
            for column in columns:
                row_list_pairs.append(cds.get_codon_pair_freq(row, column))
                row_list_expected.append(cds.get_codon_freq(row)[1]*cds.get_codon_freq(column)[1])
            dict_vals_pairs[row] = row_list_pairs
            dict_vals_expected[row] = row_list_expected
        pair_values = pd.DataFrame.from_dict(dict_vals_pairs, orient='index', columns=columns )
        expected_values = pd.DataFrame.from_dict(dict_vals_expected, orient='index', columns=columns)
        diff = np.subtract(pair_values,expected_values)

        return_file[diff_name] = diff
        sheet_names[diff_name] = diff_name
        return_file[pair_values_name] = pair_values
        sheet_names[pair_values_name] = pair_values_name
        return_file[expected_values_name] = expected_values
        sheet_names[expected_values_name] = expected_values_name

        for i in cds.CDS_info:
            dict_vals_pairs = {}
            dict_vals_expected = {}
            pair_values_name = "pair values " + i["label"].split()[0]
            expected_values_name = "expected values " + i["label"].split()[0]
            diff_name = "diff " + i["label"].split()[0]
            for row in rows:
                row_list_pairs = []
                row_list_expected = []
                for column in columns:
                    row_list_pairs.append(cds.get_codon_pair_freq(row, column, i))
                    row_list_expected.append(cds.get_codon_freq(row,i)[1]*cds.get_codon_freq(column,i)[1])
                dict_vals_pairs[row] = row_list_pairs
                dict_vals_expected[row] = row_list_expected
            pair_values = pd.DataFrame.from_dict(dict_vals_pairs, orient='index', columns=columns)
            expected_values = pd.DataFrame.from_dict(dict_vals_expected, orient='index', columns=columns)

            return_file[diff_name] = diff
            sheet_names[diff_name] = diff_name
            return_file[pair_values_name] = pair_values
            sheet_names[pair_values_name] = pair_values_name
            return_file[expected_values_name] = expected_values
            sheet_names[expected_values_name] = expected_values_name

        with pd.ExcelWriter(excel_out_name) as writer:
            for i in sheet_names:
                return_file[i].to_excel(writer, sheet_name=i)
            # pair_values.to_excel(writer, sheet_name="pair values")
            # expected_values.to_excel(writer, sheet_name="expected values")
    print("yay! finished making file!")
    return return_file

def plotting(in_gb_file, cds):
    # makes an excel file from the codon data, and the pairing of the codons
    # if an excel file with the correct name already exists, it justs retrieves that one.
    codon_pair_data_set = make_excel_file(in_gb_file, cds)
    pv = codon_pair_data_set["pair values"]
    ev = codon_pair_data_set["expected values"]

    dm = np.subtract(pv,ev)
    # nm = dm/abs(dm)
    # sq = np.power(dm, 3)
    # # ax = sns.heatmap(dm, linewidths=0.5)
    plt.imshow(dm, "seismic", vmax=0.01, vmin=-0.01)
    plt.colorbar()
    plt.savefig('figures/codon_pair.png')
    plt.show()


if __name__ == '__main__':
    # Input files
    # in_gb_file = "test-sequence.gb"
    in_gb_file = "genbank\cp.gb"
    in_fasta_file = "fasta\genes.fasta"

    # creates the object containing all the CDS infoomation in the input file
    cds = CDS(in_gb_file, in_fasta_file)
    plotting(in_gb_file, cds)

    """
    # makes an excel file from the codon data, and the pairing of the codons
    # if an excel file with the correct name already exists, it justs retrieves that one.
    codon_pair_data_set = make_excel_file(in_gb_file, cds)
    pv = codon_pair_data_set["pair values"]
    ev = codon_pair_data_set["expected values"]


    dm = np.subtract(pv,ev)
    # nm = dm/abs(dm)
    # sq = np.power(dm, 3)
    # # ax = sns.heatmap(dm, linewidths=0.5)
    plt.imshow(dm, "seismic", vmax=0.01, vmin=-0.01)
    plt.colorbar()
    plt.savefig('codon_pair.png')
    plt.show()
    """





