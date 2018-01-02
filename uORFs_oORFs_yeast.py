import re
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pickle
from intermine.webservice import Service
import kozakScoring as k
import random
import re
import score_uORFs_functions as f
import pandas as pd
import seaborn as sns

sns.set_style("dark")
print 'running....'

with open('pDL1728') as fp1:
    fp1 = fp1.readlines()

firefly_seq = fp1[0]


def Ring():
    os.system('say "your script is done"')
    
def Broken():
    os.system('say "your program is broken"')




'''
with open('optimal_codons.txt') as fp:
    fp = fp.readlines()[14:]
codons, fraction =  [], []
for i in fp:
    i = i.strip('\n')
    linespace = i.split('    ')
    if len(linespace)>1:
        codon = linespace[1].replace(" ", "")
        codons.append(codon)
        fraction.append(float(linespace[4]))

optimal_codons = zip(codons, fraction)



def ScoreCodons (seq):
    newSeq = [seq[i:i+3] for i in range(3, len(seq)-3, 3)]
    codonScore = []
    for i in newSeq:
        for optimalcodon in optimal_codons:
            if i == optimalcodon[0]:
                codonScore.append(optimalcodon[1])
    d = 10000 *(reduce(lambda x, y: x*y, codonScore))
    print d



ScoreCodons('ATGGGGTCTATCGATCGATTCGTATAA')

'''

################################################################################################
################################################################################################
######################      OPEN ALL FILES           ###########################################
################################################################################################
################################################################################################

def MakeList(fileName):
    '''returns a list of tuples containing (geneName, tl_length)'''
    Naga = []
    for n in fileName:
        gene = n.split('\t')[0]
        TL = n.split('\t')[1]
        TL = TL.strip('\n')
        item = (gene, TL)
        Naga.append(item)
    return Naga

promoters = open(r"yeast promoter sizes from YPA.txt")
crac_xell = pd.read_excel('UPF1_crac.xlsx', skiprows= 0, parse_cols = range(0,3,1))
nagala = open('Nagalakshmi_tl_lengths', 'r')
Naga = MakeList(nagala)
Air = open('Airebere_Longest_TLs', 'r')
air = f.MakeList(Air)
yeastGeneNames = open('resultsnames.tsv', 'r')
trans_rate = open(r"translation_rate_arava_et_al.txt")
trans_rate = list(trans_rate)

################################################################################################
################################################################################################
######################      Define Class           ###########################################
################################################################################################
################################################################################################

differentGeneNames = {}

for i in yeastGeneNames:
    linestring = i.rstrip()
    geneName = linestring.split('\t')[0]
    shortName = linestring.split('\t')[1]
    differentGeneNames[shortName] = geneName

class Gene:

    def __init__ (self,SystematicName, five_prime_seq, ORFseq, symbol):
        self.five_prime_seq = five_prime_seq
        self.ORFseq = ORFseq
        self.symbol = symbol
        self.SystematicName = SystematicName
        

    class uAUG:
        def __init__ (self):
            pass

    def TLseqNaga(self):
        try:
            return self.five_prime_seq[-int(self.TL_length_Naga):]
        except AttributeError:
            pass

    def TLseqAir(self):
        try:
            return self.five_prime_seq[-int(self.TL_length_Air):]
        except AttributeError:
             pass


    def Chose_TL(self, SPECIFIED_LENGTH  = 0):
        if SPECIFIED_LENGTH == 0:
            try:
                return self.five_prime_seq[-int(self.TL_length_Air):]
            except AttributeError:
                try:
                    return self.five_prime_seq[-int(self.TL_length_Naga):]
                except AttributeError:
                    pass
        else:
            return self.five_prime_seq[-int(SPECIFIED_LENGTH):]

    def TL_length_Choice (self):
        x =  self.Chose_TL()
        if x is not None:
            return len(x)
        else:
            return

    def MainORFkozak(self, s):
       if s == 'scramble':
           five_prime = shuffle_dna_seq(self.five_prime_seq)
           five_prime = five_prime[-4:]
           after_atg = shuffle_dna_seq(self.ORFseq)
           after_atg = after_atg [3:5]
        
       else:
           five_prime = self.five_prime_seq[-4:]
           after_atg = self.ORFseq[3:5]
       koz_seq = five_prime + after_atg
       
       if len(koz_seq) == 6:
           return ScoreKozak(koz_seq)
       else:
           return 0



    def uORFinfo(self, s):               
        
        TLseq = self.Chose_TL()
        
        if s == 'scramble':
            TLseq = shuffle_dna_seq(TLseq)

        if TLseq is not None:
            all_AUGs = find_ATG(TLseq)

            CLASS_LIST= []

            for i in all_AUGs:
                new = self.uAUG()
                new.start = i[0]
                new.frame = i[1]

                dist_to_mainORF = len(TLseq) - new.start
                kozSeq = AssembleSeqtoScoreKozakforuaAUG(TLseq, new.start)
                score = ScoreKozak(kozSeq)
                LEN_OF_uORF = Define_if_ORF(TLseq, new.start)
                LEN_OF_aORF = Define_if_aORF(dist_to_mainORF, self.ORFseq, TLseq, new.start)
               # LEN_OF_aORF = Define_if_aORF(dist_to_mainORF, firefly_seq, TLseq, new.start)

                new.uORF_length = LEN_OF_uORF
                new.aORF_LEN = LEN_OF_aORF
                new.KozakScore = score
                new.dist_to_mainORF = dist_to_mainORF
                CLASS_LIST.append(new)

            return CLASS_LIST



#################################################################################################
################################################################################################
######################      functions used by the class but not methods          ###############
################################################################################################
################################################################################################

def shuffle_dna_seq(seq):
    if seq is not None:
        x= list(seq)
        random.shuffle(x)
        return ''.join(x)
    else:
        return None


def ScoreKozak(SEQUENCE):
    'takes a sequence of 5 nucleotides and returns a Kozak score'
    if(SEQUENCE[0]=='A'):
        score_min1 = 0.439928
    if(SEQUENCE[0]=='T'):
        score_min1 = 0.21725
    if(SEQUENCE[0]=='G'):
        score_min1 = 0.164824
    if(SEQUENCE[0]=='C'):
        score_min1 = 0.178123
    if(SEQUENCE[1]=='A'):
        score_min2 = 0.397639
    if(SEQUENCE[1]=='T'):
        score_min2 = 0.240137
    if(SEQUENCE[1]=='G'):
        score_min2 = 0.143903
    if(SEQUENCE[1]=='C'):
        score_min2 = 0.21832
    if(SEQUENCE[2]=='A'):
        score_min3 = 0.576061
    if(SEQUENCE[2]=='T'):
        score_min3 = 0.140167
    if(SEQUENCE[2]=='G'):
        score_min3 = 0.189331
    if(SEQUENCE[2]=='C'):
        score_min3 = 0.094441
    if(SEQUENCE[3]=='A'):
        score_min4 = 0.42902
    if(SEQUENCE[3]=='T'):
        score_min4 = 0.224298
    if(SEQUENCE[3]=='G'):
        score_min4 = 0.148984
    if(SEQUENCE[3]=='C'):
        score_min4 = 0.197699
    if(SEQUENCE[4]=='A'):
        score_plus1 = 0.317693
    if(SEQUENCE[4]=='T'):
        score_plus1 = 0.261805
    if(SEQUENCE[4]=='G'):
        score_plus1 = 0.288255
    if(SEQUENCE[4]=='C'):
        score_plus1 = 0.132247
    if(SEQUENCE[5]=='A'):
        score_plus2 = 0.265093
    if(SEQUENCE[5]=='T'):
        score_plus2 = 0.216378
    if(SEQUENCE[5]=='G'):
        score_plus2 = 0.151973
    if(SEQUENCE[5]=='C'):
        score_plus2 = 0.366557
    KozakScore= (score_min1*score_min2*score_min3*score_min4*score_plus1*score_plus2)*10000

    return KozakScore


def AssembleSeqtoScoreKozakforuaAUG(seq, start):
    seq = seq+'ATG'
    ''' takes TL and the start and returns a sequence to be used in ScoreKozak'''
    seq = list(seq)
    KozakSeq = seq[start-4] + seq[start-3] + seq[start-2]  + seq[start-1] + seq[start+3] + seq[start+4]
    return KozakSeq


def find_STOP(sequence):
    STOP = []
    for m in re.finditer('TGA', sequence):
        STOP.append(m.start())
    for m in re.finditer('TAG', sequence):
        STOP.append(m.start())
    for m in re.finditer('TAA', sequence):
        STOP.append(m.start())
    STOP.sort(key=int)
    return STOP

def Define_if_ORF(seq, start):
    '''Take a number (ATG location), and a sequence and return a length if uORF. if not an ORF return 0'''

    stops = find_STOP(seq)
    if len(stops) == 0 :
        return 0
    else:
        for y in stops:
            if y > start and (y - start) % 3 == 0:
                length_of_uORF = y - start + 3
  
                break
            else:
                length_of_uORF = 0
                

    return length_of_uORF


def Define_if_aORF(bpsTomainATG, orfSeq, seq, start):
    '''take a number which the distance from the uAUG to the main ATG, and the main ORF sequence and calculate aORF length, return 0 if not an aORF'''

    uORFlEN = Define_if_ORF(seq, start)
    if uORFlEN == 0:
        stops = find_STOP(orfSeq)
        if len(stops) == 0:
            return 0
        for i in stops:
            potential_ORF = i + bpsTomainATG
            if potential_ORF % 3 == 0:
                aORF_LEN = potential_ORF + 3
                break
            else:
                aORF_LEN = 0
        return aORF_LEN
    else:
        return 0


def find_ATG(sequence):
    '''returns a list of ATG start locations (from the 5'end) and in Frame with the 3'end '''

    all, frames = [], []
    for m in re.finditer('ATG', sequence):
        all.append(m.start())
        if (len(sequence) - m.start()) % 3 == 0:
            frames.append('IN')
        else:
            frames.append('OUT')
    INFO = zip(all, frames)

    return INFO


def RerurnContentsOfClassItem(Jumbled_class):
    attrs = vars(Jumbled_class)
    return ', '.join("%s: %s" % item for item in attrs.items())

Class_dict = {}





################################################################################################
################################################################################################
######################      get info from SGD          #########################################
################################################################################################
################################################################################################
'''
service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view("secondaryIdentifier", "symbol", "length", "flankingRegions.direction",
    "flankingRegions.sequence.length", "flankingRegions.sequence.residues")

# Uncomment and edit the line below (the default) to select a custom sort order:query.add_sort_order("Gene.secondaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("Gene", "IN", "ALL_Verified_Uncharacterized_Dubious_ORFs", code = "B")
query.add_constraint("flankingRegions.direction", "=", "both", code = "C")
query.add_constraint("flankingRegions.distance", "=", "1.0kb", code = "A")
query.add_constraint("flankingRegions.includeGene", "=", "true", code = "D")


for row in query.rows():
    NAME = row["secondaryIdentifier"]
    SystematicName = row["secondaryIdentifier"]
    shortname = row["symbol"]
    seq = row["flankingRegions.sequence.residues"]
    upstreamSeq = seq[0:1000]
    downstreamSeq = seq[-1000:]
    interseq = seq[1000:]
    ORFseq = interseq[:-1000]

    NAME = Gene(SystematicName, upstreamSeq, ORFseq, shortname)
    Class_dict[NAME.SystematicName] = NAME



with open('SGDinfo.pickle', 'wb') as handle:
  pickle.dump(Class_dict, handle)

'''

with open('SGDinfo.pickle', 'rb') as handle:
    Class_dict = pickle.load(handle)




################################################################################################
################################################################################################
######################      get info from other files         ##################################
################################################################################################
################################################################################################

'''
for i in Naga:
    try:
        Class_dict[i[0]].TL_length_Naga = i[1] #[i[0]] is the systematic gene name and  [i[1]] is the TL_length
    except KeyError:
        pass

for i in air:
    try:
        Class_dict[i[0]].TL_length_Air = i[1] #[i[0]] is the systematic gene name and  [i[1]] is the TL_length
    except KeyError:
        pass


for line in promoters:
    linestring = line.rstrip()
    linelist = linestring.split('\t')
    genename = linelist[0].split (' ')[0]

    promoter_len = int(linelist[3])
    try:
        Class_dict[genename].Promoter_length = promoter_len
    except KeyError:
        pass


for i in trans_rate:
    string=i.rstrip()
    stringsplit=string.split('\t')
    genename2=stringsplit[0]
    trans = float(stringsplit[1])
    try:
        Class_dict[genename2].translation_rate = trans
    except KeyError:
   
        pass

### since some of the crac data have gene names which have multiple
# names simply trying to match up short names doesnt work so this is required
crac_xell2 = {}

for index, row in crac_xell.iterrows():
    if row[0] in differentGeneNames:
        name =  differentGeneNames[row[0]]
        crac_xell2[name] = row[1], row[2]
    else:
        crac_xell2[row[0]] = row[1], row[2]



for index, row in crac_xell2.iteritems():

    for key, item in Class_dict.iteritems():

        if index == key or index == item.symbol:
            Class_dict[item.SystematicName].crac1 = row[0]
            Class_dict[item.SystematicName].crac2 = row[1]
        else:
            pass





with open('ClassPickle.pickle', 'wb') as handle:
  pickle.dump(Class_dict, handle)

'''


with open('ClassPickle.pickle', 'rb') as handle:
    Class_dict = pickle.load(handle)



gene_name, gene_name2, Kozak_Score,  frame, length_uORF, Crac1, Crac2, TL_Length, mainKozak, distanceToMainORF, aORF_LENs  = [], [], [], [], [], [], [], [], [], [], []



################################################################################################
################################################################################################
######################      Read dict and make DFS           ###################################
################################################################################################
################################################################################################





def Percent_of_uORFs(scr):

    colNames = 'SystematicName', 'GeneName', 'Kozak_Score',  'frame', 'length_uORF', 'Crac1', 'Crac2', 'TL_Length', 'mainKozak', 'distanceToMainORF', 'aORF_LENs', 'TranslationRate'

    df11 = pd.DataFrame(columns= colNames)
    num=1
    for i, items in Class_dict.iteritems():
        num = 1+num
        miniDict = {}

        info = items.uORFinfo(scr)
        
        if info is not None and len(info) > 0:

            for x in info:

                if x.uORF_length > 0 and x.aORF_LEN > 0:
                    print 'WARNING!!!!!!!!!!!!'

                Kozak_Score = x.KozakScore
                aORF_LEN = x.aORF_LEN
                frame = x.frame
                start = x.start
                aORF_LEN = x.aORF_LEN
                length_uORF = x.uORF_length
               # if length_uORF == 6:
                #    print length_uORF, items.symbol
                distanceToMainORF = x.dist_to_mainORF
                
                GeneName = items.symbol
                mainKozak = items.MainORFkozak(scr)
                try:
                    TranslationRate = items.translation_rate
                except:
                    TranslationRate = 'n/a'

                SystematicName = i



                if items.TL_length_Choice() is not None:
                    TL_Length = items.TL_length_Choice()
                else:
                    TL_Length ='unknown'

                try:
                    crac1 = items.crac1
                    crac2 = items.crac2
                except:
                    crac1 = 0
                    crac2 = 0


                miniDict[num] = SystematicName, GeneName, Kozak_Score, frame, length_uORF, crac1, crac2, TL_Length, mainKozak, distanceToMainORF, aORF_LEN, TranslationRate
                dfNew = pd.DataFrame(miniDict, index = colNames).transpose()
               
                df11 = df11.append(dfNew, ignore_index =True)


        else:
                Kozak_Score = 0
                aORF_LEN = 0
                frame = 0
                start = 0
                length_uORF = 0
                distanceToMainORF = 0
                aORF_LEN = 0
                GeneName = items.symbol
                mainKozak = items.MainORFkozak(scr)
                try:
                    TranslationRate = items.translation_rate
                except:
                    TranslationRate = 'n/a'

                SystematicName = i


                if items.TL_length_Choice() is not None:
                    TL_Length = items.TL_length_Choice()
                else:
                    TL_Length ='unknown'


                try:
                    crac1 = items.crac1
                    crac2 = items.crac2
                except:
                    crac1 = 0
                    crac2 = 0

                miniDict[num] = SystematicName, GeneName, Kozak_Score, frame, length_uORF, crac1, crac2, TL_Length, mainKozak, distanceToMainORF, aORF_LEN, TranslationRate
                dfNew = pd.DataFrame(miniDict, index = colNames).transpose()
                df11 = df11.append(dfNew, ignore_index=True)


    df11['into_Main_ORF'] = df11['aORF_LENs'] - df11['distanceToMainORF']
    df11['uORF_stop_to_mATG'] = df11['distanceToMainORF'] - df11['length_uORF'] 
    
 #   df11.to_csv('uORF analysis results.csv')


    def HowManyGenes(bigDF, queries, queries2, dfReturn = 0):
        '''how many genes of Query2 are in Query1'''
        if queries == 0 or queries == '' or queries ==None:
            allthequeriesDF = bigDF
        else:
            allthequeriesDF = bigDF.query(queries)
        queries2DF = allthequeriesDF.query(queries2)
        # NOW WE NEED TO SEE HOW MANY UNIQUE GENES THERE ARE
        TOTAL = len(allthequeriesDF.SystematicName.unique())
        QUERIES_2 = len(queries2DF.SystematicName.unique())
        try:
            perc = (float(QUERIES_2)/ float(TOTAL))*100.0
        except ZeroDivisionError:
            perc = 0
        if dfReturn == 0 :
            return TOTAL, QUERIES_2, perc
        else:
            return TOTAL, QUERIES_2, perc, queries2DF



    def freqGenesBindUPF1(df, d):
        aORF =  HowManyGenes(df, 'aORF_LENs > 0 and into_Main_ORF == %d and TL_Length != "unknown"'%(d), 'Crac1 > 500')
        return aORF, d



    total_df = df11


    total_df['uORF'] = total_df['length_uORF'] > 0
    total_df['aORF'] = total_df['aORF_LENs'] > 0

    uORFsOnly = total_df[total_df['uORF'] == True]
    aORFsOnly = total_df[total_df['aORF'] == True]



    total_df['Koz_Main - Koz_uORF']  = total_df['mainKozak'] - total_df['Kozak_Score']


    Kozak_differences = total_df['Koz_Main - Koz_uORF'].tolist()
    TranslationRates  = total_df.filter(['length_uORF','TranslationRate'], axis=1)
    TranslationRates = TranslationRates[TranslationRates['TranslationRate'] != 'n/a']
    TranslationRates = TranslationRates[TranslationRates['length_uORF'] != 0]


    main_ORFs = total_df.filter(['SystematicName','mainKozak'], axis=1)   
    main_ORFs = main_ORFs.drop_duplicates(subset= ['SystematicName'] )
    ORF_Kozaks = main_ORFs['mainKozak'].tolist()



    uORF_lengths = total_df['length_uORF'].tolist()
    uORF_lengths = filter(lambda a: a != 0, uORF_lengths)

    aORF_lengths = total_df[total_df['aORF'] == True]
    #aORF_lengths = TL_known[TL_known['aORF'] == True]
    aORF_lengths = aORF_lengths['aORF_LENs'].tolist()


    #aORF_lengths = filter(lambda a: a != 0, aORF_lengths)

    
    uORF_dist = uORFsOnly['uORF_stop_to_mATG'].tolist()
    uORF_dist = filter(lambda a: a != 0, uORF_dist)

    aORF_atg_to_atg = aORFsOnly['distanceToMainORF'].tolist()
    aORF_atg_to_atg = filter(lambda a: a != 0, aORF_atg_to_atg)


    aORF_dist = aORFsOnly['into_Main_ORF'].tolist()
    aORF_dist = filter(lambda a: a != 0, aORF_dist)

    
    ################################################################################################
    ################################################################################################
    ######################       plotting and analysis           ###################################
    ################################################################################################
    ################################################################################################

    
    Kozak_aORFsOnly = aORFsOnly['Kozak_Score'].tolist()
    Kozak_uORFsOnly = uORFsOnly['Kozak_Score'].tolist()

    

    def GetFreqTable(df, item):
        '''return df with counts, per gene, either a uORF or an aORF'''
        uORF_freq = df.groupby(["SystematicName", item]).size()
        uORF_freq = uORF_freq.to_frame() # convert to DF
        uORF_freq.reset_index(level=0, inplace=True)
        uORF_freq.reset_index(level=0, inplace=True)
        uORF_freq = uORF_freq.rename(columns = {0L: 'Freq'})
        uORF_freq = uORF_freq[uORF_freq[item] == True]
        return uORF_freq

   
    unique_genes_with_TLs = total_df[total_df['TL_Length'] != 'unknown']
    

    unique_genes_with_TLs_total = unique_genes_with_TLs.drop_duplicates(subset= ['SystematicName'] )

    unique_genes_with_TLs = unique_genes_with_TLs.drop_duplicates(subset= ['SystematicName'] )

    average_tl_length = unique_genes_with_TLs["TL_Length"].median()
    print 'average_tl_length:', average_tl_length
    print 'len(unique_genes_with_TLs):', len(unique_genes_with_TLs)
    
    unique_genes_with_TLs = unique_genes_with_TLs.filter(['SystematicName','mainKozak'], axis=1)
    unique_genes_with_TLs = unique_genes_with_TLs.drop_duplicates(subset= ['SystematicName','mainKozak'] )



    unique_genes_with_TLs_crac = unique_genes_with_TLs_total.filter(['SystematicName','Crac1', 'Crac2', 'uORF', 'aORF'], axis=1)


    crac_total = unique_genes_with_TLs_crac['Crac1'].tolist()


    #crac_aORF_pos = unique_genes_with_TLs_crac[(unique_genes_with_TLs_crac['aORF'] == True) |  (unique_genes_with_TLs_crac['uORF'] == True)]
    crac_aORF_pos = unique_genes_with_TLs_crac[unique_genes_with_TLs_crac['aORF'] == True]
    crac_aORF_pos = crac_aORF_pos['Crac1'].tolist()

    #crac_aORF_neg = unique_genes_with_TLs_crac[(unique_genes_with_TLs_crac['aORF'] == False) & (unique_genes_with_TLs_crac['uORF'] == False)]
    crac_aORF_neg = unique_genes_with_TLs_crac[unique_genes_with_TLs_crac['aORF'] == False]
    crac_aORF_neg = crac_aORF_neg['Crac1'].tolist()


    crac_uORF_pos = unique_genes_with_TLs_crac[unique_genes_with_TLs_crac['uORF'] == True]
    crac_uORF_pos = crac_uORF_pos['Crac1'].tolist()



    crac_uORF_neg = unique_genes_with_TLs_crac[unique_genes_with_TLs_crac['uORF'] == False]
    crac_uORF_neg = crac_uORF_neg['Crac1'].tolist()




    uORF_freq = GetFreqTable(total_df, 'uORF')
    aORF_freq = GetFreqTable(total_df, 'aORF')

    num_aORFs = len(aORF_freq)
    num_uORFs = len(uORF_freq)
    num_genes_with_TL = len(unique_genes_with_TLs)

    num_uORFs, num_aORFs, num_genes_with_TL = float(num_uORFs),float(num_aORFs), float(num_genes_with_TL)

    frac_uORFs = num_uORFs/num_genes_with_TL
    frac_aORFs = num_aORFs/num_genes_with_TL



    return frac_uORFs, frac_aORFs, Kozak_uORFsOnly, Kozak_aORFsOnly, \
           ORF_Kozaks, uORF_lengths, aORF_lengths, uORF_dist, aORF_dist, total_df, \
           aORF_atg_to_atg, crac_total, crac_aORF_pos, crac_aORF_neg,\
           num_genes_with_TL, Kozak_differences




uORFs_and_aORFs_Percents = []

Kozak_uORFsOnly, Kozak_aORFsOnly, ORF_Kozaks = [], [], []
Kozak_uORFsOnly_scr, Kozak_aORFsOnly_scr, ORF_Kozaks_scr = [], [], []

#info = Percent_of_uORFs('scrambley')


for i in range(0, 1):
    info = Percent_of_uORFs('scrambleyjjjj')
 
    
    info_scramble = Percent_of_uORFs('scramble')

    
    uORF, aORF = info[0], info[1]
    uORF_scr, aORF_scr = info_scramble[0], info_scramble[1]
    print 'uORF', 'oORF'
    print uORF, aORF
 #   uORFs_and_aORFs_Percents.append(info)
    print 'uORF_scramble', 'oORF_scramble'
    print uORF_scr, aORF_scr
 #   Kozak_uORFsOnly.append(info[2])
 #   Kozak_aORFsOnly.append(info[3])
 #   ORF_Kozaks.append(info[4])



#with open('freqs', 'wb') as freqs:
#    pickle.dump(uORFs_and_aORFs_Percents, freqs)

#with open('freqs', 'rb') as freqs:
#    freq_list = pickle.load(freqs)



#info_SCR = Percent_of_uORFs('scramble')
#Kozak_differences = info[15]
#Kozak_differences_scr = info_SCR[15]


#with open('Kozak_differencesl', 'wb') as f:
 #   pickle.dump(Kozak_differences, f)

#with open('Kozak_differences_scr', 'wb') as f:
 #   pickle.dump(Kozak_differences_scr, f)



#uORF_lengths = info[5], info[14]
#aORF_lengths = info[6]




#uORF_lengths_scr = info_SCR[5], info[14]
#aORF_lengths_scr = info_SCR[6]

uORF_dist = info[7]
aORF_dist = info_scramble[7]


with open('uORF_dist', 'wb') as f:
    pickle.dump(uORF_dist, f)

with open('aORF_dist', 'wb') as f:
    pickle.dump(aORF_dist, f)




'''
aORF_atg_to_atg = info[10]
aORF_atg_to_atg_scr = info_SCR[10]



crac_total, crac_aORF_pos, crac_aORF_neg = info[11], info[12], info[13]
print len(crac_total), len(crac_aORF_pos), len(crac_aORF_neg )

with open('crac_total', 'wb') as f:
    pickle.dump(crac_total, f)

with open('crac_aORF_pos', 'wb') as f:
    pickle.dump(crac_aORF_pos, f)

with open('crac_aORF_neg', 'wb') as f:
    pickle.dump(crac_aORF_neg, f)

with open('aORF_atg_to_atg', 'wb') as f:
    pickle.dump(aORF_atg_to_atg, f)

with open('aORF_atg_to_atg_scr', 'wb') as f:
    pickle.dump(aORF_atg_to_atg_scr, f)

with open('uORF_lengths_scr', 'wb') as f:
    pickle.dump(uORF_lengths_scr, f)

with open('aORF_lengths_scr', 'wb') as f:
    pickle.dump(aORF_lengths_scr, f)


with open('uORF_lengths', 'wb') as f:
    pickle.dump(uORF_lengths, f)

with open('aORF_lengths', 'wb') as f:
    pickle.dump(aORF_lengths, f)




with open('uORF_dist_scr', 'wb') as f:
    pickle.dump(uORF_dist_scr, f)

with open('aORF_dist_scr', 'wb') as f:
    pickle.dump(aORF_dist_scr, f)


with open('uORF_dist', 'wb') as f:
    pickle.dump(uORF_dist, f)

with open('aORF_dist', 'wb') as f:
    pickle.dump(aORF_dist, f)



'''''

'''

for i in range(0, 1):
    info = Percent_of_uORFs('scrambley')
    info_scramble = Percent_of_uORFs('scramble')
    uORF, aORF = info[0], info[1]
#   print uORF, aORF
    uORFs_and_aORFs_Percents.append(info)
    Kozak_uORFsOnly.append(info[2])
    Kozak_aORFsOnly.append(info[3])
    ORF_Kozaks.append(info[4])
    


Kozak_uORFsOnly_scr = Percent_of_uORFs('scramble')[2]
Kozak_aORFsOnly_scr= Percent_of_uORFs('scramble')[3]

ORF_Kozaks_scr = Percent_of_uORFs('scramble')[4]


ORF_Kozaks = Percent_of_uORFs('c')[4]

Kozaks_main = zip( ORF_Kozaks_scr, ORF_Kozaks)

print len(ORF_Kozaks_scr), len(ORF_Kozaks), 'ooo'



with open('Kozak_uORFsOnly', 'wb') as f:
    pickle.dump(Kozak_uORFsOnly, f)


with open('Kozak_aORFsOnly', 'wb') as f:
    pickle.dump(Kozak_aORFsOnly, f)


with open('ORF_Kozaks', 'wb') as f:
    pickle.dump(ORF_Kozaks, f)
    
with open('Kozak_uORFsOnly', 'rb') as f:
    Kozak_uORFsOnly = pickle.load(f)

with open('Kozak_aORFsOnly', 'rb') as f:
    Kozak_aORFsOnly = pickle.load(f)

with open('ORF_Kozaks', 'rb') as f:
    ORF_Kozaks = pickle.load(f)



with open('ORF_Kozaks_scr', 'wb') as f:
    pickle.dump(ORF_Kozaks_scr, f)
    
'''


#with open('freqs', 'wb') as freqs:
#    pickle.dump(uORFs_and_aORFs_Percents, freqs)

#with open('freqs', 'rb') as freqs:
#    freq_list = pickle.load(freqs)



#Ring()


#writer = pd.ExcelWriter('main_ORFs.xlsx')
#n.to_excel(writer,'Sheet3')
#writer.save()


#writer = pd.ExcelWriter('output_df.xlsx')
#df11.to_excel(writer,'Sheet1')
#writer.save()
