#################################
#    Kristoffer Holm Hundvin    #
#    Master thesis              #
#################################

from dna_features_viewer import GraphicFeature, GraphicRecord
import pandas as pd

import Bio.SeqUtils
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm
import PIL.Image

from graphics import *


# draws the full gene ma of all the CDS in the genbank file
# using the GenomeDiagram
def draw_full_gene_map(gb):
    gd_diagram = GenomeDiagram.Diagram("Chlamydomonals genes")
    scale_track = gd_diagram.new_track(1, name="scale", scale_smalltick_interval=20000,
                                       scale_smalltick_labels=1, scale_fontsize=8,
                                       scale_color=colors.black, start=0, end=len(gb))

    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features", scale_smalltick_interval=20000,
                                                 scale_smalltick_labels=1, scale_fontsize=8,
                                                 scale_color=colors.black, height=2)
    gd_feature_set = gd_track_for_features.new_set()
    for feature in gb.features:
        if feature.type != "CDS":
            # Exclude this feature
            continue
        # this %2 makes them alternate between blue and light blue
        # if len(gd_feature_set) % 2 == 0:
        # using feature strand, makes the forward genes blue and rev. light blue
        if feature.strand == 1:
            color = colors.lightgreen
        else:
            color = colors.seagreen
        gd_feature_set.add_feature(feature, color=color, label=False, sigil="BIGARROW", arrowshaft_height=1.0,
                                   arrowhead_length=0.2)
        # gd_feature_set.add_feature(feature, sigil="ARROW", color="brown",
        #                            arrowshaft_height=1.0)

        gd_diagram.draw(format="linear", orientation="landscape", pagesize=(4 * cm, 29 * cm),
                        fragments=1, start=0, end=len(gb), track_size=0.05, tracklines=0, )

        gd_diagram.write("plasmid_linear.png", "PNG")


def draw_genes_w_DFV(gb_infile, full_genome=True):
    owd = os.getcwd()
    gfeatures = []
    codon_features = []
    firstCDS = True
    gb = SeqIO.read(gb_infile, "genbank")
    image_names = []
    for feature in gb.features:
        if feature.type == "CDS":
            if firstCDS:
                name = feature.qualifiers["label"][0]
                firstCDS = False
            else:
                name = feature.qualifiers["label"][0]
            start = feature.location.start.position
            end = feature.location.end.position
            sense = feature.strand
            gfeature = GraphicFeature(start=start, end=end, strand=sense, label=name, color="#cffccc")
            # gfeature.color = "#cffccc"
            gfeatures.append(gfeature)
        elif feature.type == "modified_base":
            codon_features.append(feature)
    firstCDS = True
    if(len(gfeatures) >= 5):
        for i in gfeatures:
            if(firstCDS):
                firstCDS = False
            else:
                i.label = None
    gfeatures[-1].label = name
    record = GraphicRecord(sequence_length=len(gb), features=gfeatures)
    ax, _ = record.plot(figure_width=10)
    if (len(gfeatures)>= 5):
        ax.figure.text(0.005, 0.2, "Chloroplast")
    filename = gb_infile.replace("genbank/", "")
    filename = filename.replace("extractions/", "")
    path = os.path.join("figures", "figure_" + filename.replace(".gb", ""))
    if (not os.path.exists(path)):
        os.makedirs(path)
    os.chdir(path)
    target_file = filename.replace(".gb", ".png")
    image_names.append(target_file)
    ax.figure.savefig(target_file)

    codon_labels = []
    for feature in codon_features:
        label = feature.qualifiers["label"][0]
        if label not in codon_labels:
            codon_labels.append(label)
    counter = 0
    for label in codon_labels:
        gfeatures_codons = []
        for feature in codon_features:
            current_label = feature.qualifiers["label"][0]
            if current_label == label:
                name = None
                start = feature.location.start.position
                end = feature.location.end.position
                gfeature = GraphicFeature(start=start, end=end, strand=0, label=name)
                gfeature.linewidth = 0.2
                gfeature.linecolor = "#db393c"
                # gfeature.color = "#cffccc"
                gfeatures_codons.append(gfeature)
        record = GraphicRecord(sequence_length=len(gb), features=gfeatures_codons)
        ax, _ = record.plot(figure_width=10, with_ruler=False, figure_height=0.5)
        ax.figure.text(0.005, 0.1, str(len(gfeatures_codons)) + " " + label)
        target_file = label.replace(" -> ", "to") + "_" + filename.replace(".gb", ".png")
        image_names.append(target_file)
        ax.figure.savefig(target_file)
        counter = counter + 1

    imgs = [PIL.Image.open(i) for i in image_names]

    min_img_width = min(i.width for i in imgs)
    total_height = 0
    for i, img in enumerate(imgs):
        # If the image is larger than the minimum width, resize it
        if img.width > min_img_width:
            imgs[i] = img.resize((min_img_width, int(img.height / img.width * min_img_width)), Image.ANTIALIAS)
        total_height += imgs[i].height

    img_merge = PIL.Image.new(imgs[0].mode, (min_img_width, total_height))
    y = 0
    for img in imgs:
        img_merge.paste(img, (0, y))

        y += img.height
    img_merge.save('combined_figures.png')

    os.chdir(owd)


def make_figures(infile):
    draw_genes_w_DFV(infile)


def aa_three_to_one(three_letter_amino_acid):
    return Bio.SeqUtils.IUPACData.protein_letters_3to1[three_letter_amino_acid.title()]


def aa_one_to_three(one_letter_amino_acid):
    return Bio.SeqUtils.IUPACData.protein_letters_1to3[one_letter_amino_acid.upper()]


def aa_full_name(aa):
    if (len(aa) == 1):
        aa = aa_one_to_three(aa)
    amino_acids = {'ala': 'alanine', 'arg': 'arginine', 'asn': 'asparagine', 'asp': 'aspartic acid',
                   'asx': 'asparagine', 'cys': 'cysteine', 'glu': 'glutamic acid', 'gln': 'glutamine',
                   'glx': 'glutamine', 'gly': 'glycine', 'his': 'histidine', 'ile': 'isoleucine',
                   'leu': 'leucine', 'lys': 'lysine', 'met': 'methionine', 'phe': 'phenylalanine',
                   'pro': 'proline', 'ser': 'serine', 'thr': 'threonine', 'trp': 'tryptophan',
                   'tyr': 'tyrosine', 'val': 'valine', '*' : "stop", 'stop' : "stop", 'end' : "stop"}
    if (aa.lower() in amino_acids):
        return (amino_acids[aa.lower()].title())
    else:
        print("unknown amino acid")
        return ("unknown amino acid")

def get_base_change_pos(c1,c2):
    changes = []
    if(len(c1)==0 or len(c2)==0):
        return []
    for i in range(len(c1)):
        if(c1[i] != c2[i]):
            changes.append(i)
    return changes


def make_codon_box(name, startPoint, codon_list, sizes):



    num_smallBox = len(codon_list)
    figures = []
    figdict = {}
    # smallBoxHeight = 20
    # smallBoxVerticalSpacing = 10
    # boxHeight = (smallBoxHeight+smallBoxVerticalSpacing) * num_smallBox + 10
    # boxWidth = win.width / 6
    # smallBoxWidth = boxWidth/3
    # smallBoxHorizontalSpacing = (boxWidth/2)-(smallBoxWidth/2)

    (smallBoxHeight,
    smallBoxVerticalSpacing,
    boxHeight,
    boxWidth,
    smallBoxWidth,
    smallBoxHorizontalSpacing) = sizes

    lightBlue = color_rgb(230,250,255)
    lightGreen = color_rgb(180, 250, 180)
    lightRed = color_rgb(250,210,215)
    lightYellow = color_rgb(255,255,210)

    #textBox = Rectangle(startPoint,Point(startPoint.x+50,startPoint.y+20))
    namePoint = Point(startPoint.x+boxWidth/2,startPoint.y+10)
    nameText = Text(namePoint,name)#textBox.getCenter(),name)
    nameText.setStyle("bold")
    nameText.setFace("courier")

    #figures.append(textBox)
    figures.append(nameText)
    startPoint = Point(startPoint.x,startPoint.y+20)


    secondPoint = Point(startPoint.x+boxWidth,startPoint.y+boxHeight)

    backBox = Rectangle(startPoint, secondPoint)
    backBox.setWidth(2)
    backBox.setFill(lightBlue)
    # backBox.draw(win)
    figures.append(backBox)
    for i in range(num_smallBox):
        #pointone = Point(startPoint.x+boxWidth/3,startPoint.y + (boxCount*smallBoxHeight)+10)
        x1 = startPoint.x+smallBoxHorizontalSpacing
        y1 = startPoint.y+((smallBoxHeight+smallBoxVerticalSpacing)*i)+smallBoxVerticalSpacing
        pointOne = Point(x1,y1)
        x2 = pointOne.x+smallBoxWidth
        y2 = pointOne.y+smallBoxHeight#+(smallBoxHeight*i)
        pointTwo = Point(x2,y2)
        rec = Rectangle(pointOne,pointTwo)
        rec.setWidth(2)


        # For a gradient when coloring
        # r = 0
        # g = 0
        # freq = codon_list[i][2]
        # if(freq >= 10):
        #     r = max((250-(freq-10)*7),180)
        #     g = 250
        # else:
        #     r =250
        #     g = min(180 +(freq*7),250)
        # lightGreen = color_rgb(int(r),int(g),180)

        # set strict valuesd for HF and LF, 10 and over are HF and 5 and lower are LF
        if codon_list[i][2] >= 10:
            rec.setFill(lightGreen)
        elif codon_list[i][2] <= 5:
            rec.setFill(lightRed)
        else:
            rec.setFill(lightYellow)
        # rec.draw(win)
        figures.append(rec)
        textPoint = Point(((x1+x2)/2), ((y1+y2)/2))



        midx = (x1+x2)/2
        middiff = (midx-x1)/2

        b1 = Point(midx-middiff, ((y1 + y2) / 2))
        b2 = Point(midx, ((y1 + y2) / 2))
        b3 = Point(midx+middiff, ((y1 + y2) / 2))
        t1 = Text(b1,codon_list[i][0][0])
        t2 = Text(b2,codon_list[i][0][1])
        t3 = Text(b3, codon_list[i][0][2])

        t1.setStyle("bold")
        t2.setStyle("bold")
        t3.setStyle("bold")

        text_lst = [t1,t2,t3]
        for j in get_base_change_pos(codon_list[i][0],codon_list[i][1]):
            text_lst[j].setTextColor("red")

        text = Text(textPoint,codon_list[i][0])
        #figures.append(text)
        figures.append(t1)
        figures.append(t2)
        figures.append(t3)

        figdict[codon_list[i][0]] = text
    arrowCount = 1
    for i in codon_list:
        if(i[1] != ""):
            arrows = make_arrow_from_box_to_box(figdict[i[0]], figdict[i[1]],
                                                arrowCount,smallBoxWidth, smallBoxHorizontalSpacing)
            arrowCount = arrowCount+1
            for i in arrows:
                figures.append(i)
    return figures

def make_arrow_from_box_to_box(c1,c2,n, textWidth,smallBoxHorizontalSpacing = 0):

    p1 = c1.getAnchor()
    p2 = c2.getAnchor()
    y1 = p1.y
    y2 = p2.y
    if(n == 1 or n == 2):
        horizontal_line_width = ((3-n)*smallBoxHorizontalSpacing)/3
        x1 = p1.x-(textWidth/2)
        x2 = x1-horizontal_line_width
    elif(n == 3 or n == 4):
        if(n == 3):n = 4
        else:n = 3
        horizontal_line_width = ((n-2)*smallBoxHorizontalSpacing/3)
        x1 = p1.x+(textWidth/2)
        x2 = x1+horizontal_line_width

    p0 = Point(x1, y1)
    p1 = Point(x2,y1)
    p2 = Point(x2,y2)
    p3 = Point(x1, y2)
    topLine = Line(p0,p1)
    line = Line(p1,p2)
    botLine = Line(p2,p3)
    botLine.setArrow("last")
    return(topLine,line,botLine)

def make_recodeing_fig(excel_file):
    fil = pd.read_excel(excel_file)
    codons = {}
    for i, j in fil.iterrows():
        AA = j[fil.columns[0]]
        codon = j[fil.columns[1]]
        newcodon = j[fil.columns[2]]
        freq = j[fil.columns[3]]
        if (pd.isnull(newcodon)):
            newcodon = ""
        if (not pd.isnull(codon)):
            if (AA in codons):
                codons[AA].append([codon.replace("*", ""), newcodon.replace("*", ""),freq])
            else:
                codons[AA] = [[codon.replace("*", ""), newcodon.replace("*", ""),freq]]

    for i in codons.keys():
        codons[i] = sorted(codons[i], key=lambda x: -len(str(x)))

    codons= {k: v for k, v in sorted(codons.items(), key=lambda item: len(item[1]))}

    winWidth = 794
    winHeight = 1123
    win = GraphWin("test", winWidth, winHeight)
    win.setBackground("white")
    figures = []
    figures_big = []

    x = 20
    y = 0
    spaceing = 20
    switch = 0

    # to start the first time on 0
    lineNr = False
    prevHeight = 0
    for i in codons:
        if (len(codons[i]) > 2):
            num_smallBox = len(codons[i])
            smallBoxHeight = 20
            smallBoxVerticalSpacing = 10
            boxHeight = (smallBoxHeight + smallBoxVerticalSpacing) * num_smallBox + 10
            boxWidth = win.width / 6
            smallBoxWidth = boxWidth / 3
            smallBoxHorizontalSpacing = (boxWidth / 2) - (smallBoxWidth / 2)

            sizes = [smallBoxHeight,smallBoxVerticalSpacing,boxHeight,boxWidth,smallBoxWidth,smallBoxHorizontalSpacing]

            if (len(codons[i])!=switch or x >= winWidth):
                switch = len(codons[i])
                y = (prevHeight + y+20)*lineNr
                lineNr = True
                x = 20
            else:
                x = x + boxWidth + spaceing
            point = Point(x,y)
            fig = make_codon_box(aa_full_name(i), point, codons[i], sizes)
            figures_big = figures_big+fig
            figures.append(fig)
            prevHeight = boxHeight
    y = y+30
    spaceing = 5
    for i in codons:
        if (len(codons[i]) <= 2):

            num_smallBox = len(codons[i])
            smallBoxHeight = 20
            smallBoxVerticalSpacing = 10
            boxHeight = (smallBoxHeight + smallBoxVerticalSpacing) * num_smallBox + 10
            boxWidth = 66.16# win.width / 6
            smallBoxWidth = 44.1#/boxWidth / 3
            smallBoxHorizontalSpacing = 10

            sizes = [smallBoxHeight,smallBoxVerticalSpacing,boxHeight,boxWidth,smallBoxWidth,smallBoxHorizontalSpacing]

            if (len(codons[i])!=switch or x >= winWidth):
                switch = len(codons[i])
                y = (prevHeight + y+20)*lineNr
                lineNr = True
                x = 20
            else:
                x = x + boxWidth + spaceing
            point = Point(x,y)
            fig = make_codon_box(i, point, codons[i], sizes)
            figures_big = figures_big+fig
            figures.append(fig)
            prevHeight = boxHeight

    for i in figures:
        for j in i:
            j.draw(win)

    win.getMouse()  # pause for click in window


    # saves the current TKinter object in postscript format
    win.postscript(file="figures/reassigned_codons.eps", colormode='color')

    win.close()


if __name__ == '__main__':
    infile = "genbank/cr_cp.gb"
    infile_cr = "genbank/cr_cp2507.gb"
    # seq_1 = load_record(os.path.join("sequences", "test-sequence2.gb"))
    # seq_2 = load_record(os.path.join("sequences", "jan_test_CR.gb"))
    #
    # diff_blocks = DiffBlocks.from_sequences(seq_1, seq_2)
    # ax1, ax2 = diff_blocks.plot(figure_width=8)
    # #ax1.figure.savefig("diff_blocks.png", bbox_inches='tight')
    # ax1.set_title("asvvdasd")
    # ax1.text(0.5,0.5,"asd")
    # ax1.figure.show()
    # gb1 = SeqIO.read(infile, "genbank")
    # gb2 = SeqIO.read(infile_cr, "genbank")

    draw_genes_w_DFV(infile)

    draw_genes_w_DFV("genbank/extractions/ex_rbcL_to_psbI.gb")
    draw_genes_w_DFV("genbank/extractions/ex_psbH_to_psbB.gb")

    draw_genes_w_DFV("genbank/extractions/ex_rbcL.gb")
    draw_genes_w_DFV("genbank/extractions/ex_atpA.gb")
    draw_genes_w_DFV("genbank/extractions/ex_psbI.gb")

    draw_genes_w_DFV("genbank/extractions/ex_psbB.gb")
    draw_genes_w_DFV("genbank/extractions/ex_psbT.gb")
    draw_genes_w_DFV("genbank/extractions/ex_psbN.gb")
    draw_genes_w_DFV("genbank/extractions/ex_psbH.gb")


    # change to true to make the reassigment box figures
    if True:
        make_recodeing_fig("excel/codon_reassignment.xlsx")

