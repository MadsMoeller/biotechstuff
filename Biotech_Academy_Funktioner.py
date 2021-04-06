# Funktioner er lavet til Biotech Academy projeketet: Regnskov + Supercomputer = Fremtidens medicin. 
# Funktionerne og materialet præsenteret er ikke akademisk korrekt, men er blevet ændret for undervisnings formål. 

# Hotfix til at kompensere for at Bio.Alpabet ikke længere understøttes.
# Det er bare statementet "from Bio.Alphabet import IUPAC" (linje 14), som er blevet udkommenteret.
# Umiddelbart ser det ud til at have virket, og at opgaverne stadigvæk kan løses.

# ==============================================================================

#import nødvendige pakker
from IPython.display import Image
from Bio import SeqIO
from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import requests
import json
from bioservices import KEGG
k = KEGG(verbose=False)
from dna_features_viewer import BiopythonTranslator
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord   
from IPython.display import YouTubeVideo
from IPython.display import HTML
from google.colab import files
from IPython.display import Image
import warnings
from Bio import BiopythonParserWarning
warnings.simplefilter('ignore', BiopythonParserWarning)

def Information_fra_GenBank_fil(information):
    record = SeqIO.read('Artemisia%20annua.gb', "genbank")
    #print(record.annotations.keys())

    seq_info = record.annotations['structured_comment']
    for k, v in seq_info.items():
      #print(v['Sequencing Technology'])
        v['Sequencing Technology'] = 'PacBio & Illumina'
        Coverage = v['Genome Coverage'] #180x er det meget eller lidt i dette tilfælde?
        sekventerings_teknologi = v['Sequencing Technology']

    if information == "organisme":
        organisme = record.annotations['organism']
        print(organisme)

    #Hvilken slags data er det?
    elif information == "molecule_type":
        print(record.annotations['molecule_type'])

    elif information == "Coverage":
        print(Coverage)

    elif information == "sekventerings_teknologi":
        print(sekventerings_teknologi)

    #Hvilken type sekventerings type er dataen lavet med?
    elif information == "sekventerings_type":
        print(record.annotations['keywords'])

    #Hvad er taxonomyen?
    elif information == "taxonomy":
        print(record.annotations['taxonomy'])

    #Print al information
    elif information == "al_information":
        print(record)

    #Print DNA
    elif information == "sekvens":
        print('Totale antal nukleotider:',len(record.seq)) #record.seq[0:10000]
        print('Antal ukendte nukleotider: 16158 \n')
        print('Eksempel på sekvens der indeholder gap:')
        print('TTGTTGATATGGAGTATTTATCCCTGTGTCNNNNNNNNNNNTACGGTTTGAAGACTCAGGAAACTCTCATTAAGCGATCAACGTAGCATGATCATCAAAAGCATGGTTTTGTAAAATCCAAGTATCCAAGCAATTGGTTCACCCTTTCAACATCCAAGTATCCAAGCAATTGGTTCACCCTTTCAACATCCAAGTATCCAAGCAATTGGTTCACCCTTTCAACATCCAAGTATCCAAGCAATTGGTTCACCCTTTCAACATCCAAGTATCCAAGCAATTGGTTCACCCTTTCAACCTCGNNNNNNNNNNNNNNNNNNN')

    else:
        print("Du har indtastet informations typen forkert. Prøv igen.")

    return

def display_yotube_video(url, **kwargs):
    id_ = url.split("=")[-1]
    return YouTubeVideo(id_, **kwargs)

def Info_om_gener(Gen_navn):
    import warnings
    from Bio import BiopythonParserWarning
    warnings.simplefilter('ignore', BiopythonParserWarning)

    file_navn = 'Artemisia%20annua.gb'

    #Sæt faste parameter
    qualifier = "gene"
    feature_type = "CDS"
    x = 0

    if Gen_navn == 'alle':
      for record in SeqIO.parse(open('Artemisia%20annua.gb',"r"), "genbank"):
        print("Antal protein kodende sekvenser:", len(record.features),"\n")
        for i,feat in enumerate (record.features):
            if (feat.type == "CDS"):
                product = feat.qualifiers['product'][0]
                print("Gen index: %d   Funktion: %s" % (i, product))

    #Find genet og print den rigtige information
    else:
        for gb_record in SeqIO.parse(open(file_navn,"r"), "genbank"):
            for (index, feature) in enumerate(gb_record.features) :
                if feature.type=="CDS":
                    if "gene" in feature.qualifiers :
                        for value in feature.qualifiers["gene"]:
                            if value == Gen_navn:
                                print("Fundet genet:",value,2*"\n","Information om genet:","\n",feature,"\n") #"Index nr:",index
                                x = 1
        if x != 1:
            print("Genet er ikke fundet")
          
    return


def sekvens_transformer_funktion(ønsket_sekvens, index):
    import warnings
    from Bio import BiopythonParserWarning
    warnings.simplefilter('ignore', BiopythonParserWarning)
    
    for gb_record in SeqIO.parse(open('Artemisia%20annua.gb',"r"), "genbank"):
        gb_feature = gb_record.features[index]
        coding_dna = gb_feature.extract(gb_record.seq) #Hent codende dna for gen
        template_dna = coding_dna.reverse_complement() #Se den modsatte streng
        messenger_rna = coding_dna.transcribe() #Konverter til mRNA
        protein_sequence = messenger_rna.translate() #Konverter til protein sekvens

        if ønsket_sekvens == "dna":
            print(coding_dna)

        elif  ønsket_sekvens == "reverse_dna":
            print(template_dna)

        elif ønsket_sekvens == "mRNA":
            print(messenger_rna)

        elif ønsket_sekvens == "protein":
            print(protein_sequence)

        else:
            print("Du har ikke indtastet typen af den ønsket sekvens korrekt. Prøv venligst igen")
    return


class ChangeFeatures(BiopythonTranslator):
  def compute_feature_color(self, feature):
      if feature.type == "CDS":
          return "pink"
      elif feature.type == "assembly_gap":
          return "blue"
      #elif feature.type == "mRNA":
      #    return "gold"
      #elif feature.type == "gene":
      #    return "green"
      else:
          return "white"

  def compute_feature_label(self, feature):
      if feature.type == 'assembly_gap':
          return "Gap"
      elif feature.type == "CDS":
          return "Gen"
      #else:
      #    return BiopythonTranslator.compute_feature_label(feature)

  def compute_filtered_features(self, features):
      #"""Vis ikke andet end CDS og gaps"""
      #qualifier = "locus_tag" #codon_start, locus_tag, protein_id
      feature_type = "CDS" #mRNA, gene, gap_assembly...

      return [feature for feature in features if (feature.type == "CDS") or (feature.type == "assembly_gap")] #or ("CT" in str(feature.qualifiers[qualifier]))]


def Visualiser_sekvens(gen):
    import warnings
    from Bio import BiopythonParserWarning
    warnings.simplefilter('ignore', BiopythonParserWarning)
    
    if gen == 'alle':
        fil = 'Artemisia%20annua.gb' 
        graphic_record = ChangeFeatures().translate_record(fil)
        ax, _ = graphic_record.plot(figure_width=20)
        ax.figure.tight_layout()

    elif gen == "aldh1":
        sequence = "CTGTGTCTAGATTTACGGTTTTGTTGAGTATGGAGTATTTATCCCTGTGTCTAGATTTACGGTTTGAAGACTCAGGAAACTCTCATTAAGCGATCAACGTAGCATGATCATCAAAAGCATGGTTTTGTAAACTCGACATGTCAATGTACCAGCCGATCCAAGTATCCAAGCAATTGGTTCACCACACCAAAAGAGTTTTACACTTAAAAACAACAATTAATTCTAAATAGTCTATGTAATGAAATATGTTTTGTGTGGGTTAGTTTAGTTCATAGTTGCGCCATAAGTATTTACAGCAA"
        record = GraphicRecord(sequence=sequence, features=[
            GraphicFeature(start=0, end=28, strand=+1, color='#ffd700',label="Promotor"),
            GraphicFeature(start=29, end=299, strand=+1, color="#ffcccc",label="aldh1")
        ])
        ax, _ = record.plot(figure_width=50)
        record.plot_sequence(ax)
        record.plot_translation(ax, (29, 299), fontdict={'weight': 'bold'})

    elif gen == 'CYP71AV1':
        sequence = "ATTTTTGGGGGCCCCCCCCCATTTTTTGGGGGGCGCGCGATGAAGTTGGTCATTCGAAATATACTTCCAAAATATGAAGTTGGTCATTCGAAATATACTTCCAAACAACCGAGCTGGTCAGGTAGATTTTGTTTCAGATGAAGATGCAATCCACCGTTGGGGGAGTTTCATGAATAACAATCGCAAATAAGATATATTGTTGATTCTTGATGATGTTTGGTCTGATACCATCATCACCGACCTCCAATTCAGGTCACGTGGATACAAGATCCTCGTGACCTCTGAAACAACCTTTAAGAGATTCGATACATATAAAGTGAGACCTCTCAGTGTTCAAGATGCCATCAATCTGTTATGCTATTCAACACTTTCGGAGCGTGCAAGTCAAGCCACAAATGACATACAGACCTTGTTGACAAGGTGAAATTTCAAATTATTCCAAGATTCATGTTTCATACCTTTATAAGAAAGTAATATCTAAACCATATTAACAAATACTAACAATTAACTTTCAAATGTTTTTGTAGTTAACCAAATGTTGCAAGAAGAATCCGCTCGCCTTAAGTGTCATTGGTGGTCGCCTAAAGGGGACACAAATGGAAAGTTGGCATCATACACTGAAAAAGCTATCTCAAGCCACACACCCTCTTATCGACCTTCCTTTGGATGAGGCAAACAGATTTCATCTCGCAAGAGCTCTCGGTTTACTCAAAGATGATGAACGCAACAGCCCCAGAAGTTCAACCTCGAAATTGACCCGATCTTACCAAGTCA"
        record = GraphicRecord(sequence=sequence, features=[
            GraphicFeature(start=1, end=38, strand=+1, color='#cffccc',label="Promotor"),
            GraphicFeature(start=39, end=774, strand=+1, color="#cff77d",label="CYP71AV1")
        ])
        ax, _ = record.plot(figure_width=100)
        record.plot_sequence(ax)
        record.plot_translation(ax, (39, 774), fontdict={'weight': 'bold'})
    return

