import score_uORFs_functions as ORF
import ORFs as o
import re
import pandas as pd
import random

#################################################################################################
################################################################################################
######################      functions used by the class but not methods START         ###############
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

    '''take a number which the distance from the uAUG to the main ATG,
    and the main ORF sequence and calculate aORF length, return 0 if not an aORF'''

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
    '''returns a list of ATG start locations (from the 3'end)  '''
    atgs = []
    for m in re.finditer('ATG', sequence):
        atgs.append(m.start())
    uATG_from_start = [-(len(sequence)-u)for u in atgs]
    return uATG_from_start


def pos_of_start(up, ORF):

    all_atgs = find_ATG(up)
    entireSeq = up+ORF
    STOPS = []
    for start in all_atgs:
        stops = find_STOP(entireSeq)
        stops = [-(len(up)-u)for u in stops]
        if len(stops) == 0 :
            STOPS.append(0)
        else:
            for stp in stops:
                if start < stp and (start - stp) % 3 == 0:
                    stopPos = stp
                    break
                else:
                    stopPos = 0
        STOPS.append(stopPos)
    if  len(all_atgs) == 0:
        all_atgs, STOPS = [0], [0]
    return zip(all_atgs, STOPS)

def RerurnContentsOfClassItem(Jumbled_class):
    attrs = vars(Jumbled_class)
    return ', '.join("%s: %s" % item for item in attrs.items())


#################################################################################################
################################################################################################
######################      functions used by the class but not methods END         ###############
################################################################################################
################################################################################################

#################################################################################################
################################################################################################
######################                 set up the classes START                   ###############
################################################################################################
################################################################################################

Class_dict = {}

class Gene:

    def __init__ (self, transcriptID, GeneID, GeneName, ORFseq):

        self.GeneID = GeneID
        self.GeneName = GeneName
        self.transcriptID = transcriptID
        self.ORFseq = ORFseq




    class uAUG:
        def __init__ (self):
            pass




    def uORFstart_stop(self, s):

        def build_uAUG(self, TLseq, ORFseq, MSG):

            start_stop = pos_of_start(TLseq, ORFseq)
            CLASS_LIST = []

            for i in start_stop:
                new = self.uAUG()
                new.start = i[0]
                new.stop = i[1]
                new.TL = MSG
                CLASS_LIST.append(new)

            return CLASS_LIST

        def build_empty_class(self, MSG):
               CLASS_LIST= []
               self.utrSeq = 'Sequence unavailable'
               new = self.uAUG()
               new.start = 0
               new.stop = 0
               new.TL = MSG
               CLASS_LIST.append(new)
               return CLASS_LIST
        try :
            TLseq = self.utrSeq
        except AttributeError:

            TLseq = 'Sequence unavailable'

        try:
            coding = self.coding
        except AttributeError:
            coding = 'Sequence unavailable'


        if TLseq  != 'Sequence unavailable' and len(TLseq) > 8 and TLseq in self.ORFseq :


           ORFseq  = self.ORFseq[len(TLseq):]
           start =  ORFseq[0:3]

           if s == 'scramble':
               TLseq = shuffle_dna_seq(TLseq)
               ORFseq = shuffle_dna_seq(ORFseq)

           if start == 'ATG':
               info = build_uAUG(self, TLseq, ORFseq, 'TL and ORF matched')
               return info

           elif start != 'ATG':
               info = build_empty_class(self, 'start not atg')
               self.utrSeq = 'Sequence unavailable'
               return info
           else:
               print 'another issue'

        elif coding in self.ORFseq and len(self.ORFseq) > len(self.coding) :


            t = re.finditer(coding, self.ORFseq)

            for it in t:
               start2, end = it.start(), it.end()

               if start2 >= 5:
                   ORFseq = self.ORFseq[start2:end]
                   TLseq_by_ORF = self.ORFseq [:start2]

                   if s == 'scramble':
                       ORFseq = shuffle_dna_seq(ORFseq)
                       TLseq_by_ORF = shuffle_dna_seq(TLseq_by_ORF)
                   info = build_uAUG(self, TLseq_by_ORF, ORFseq, 'coding seq and cDNA matched')
                   self.utrSeq = TLseq_by_ORF
                   return info

               elif start2 < 5:
                   build_empty_class(self, 'start less than 5')
                   self.utrSeq = 'Sequence unavailable'

               else:
                   print 'problem', start2


        elif TLseq  == 'Sequence unavailable':
            info = build_empty_class(self, 'TL seq listed as unavailable')
            return info

        elif TLseq not in self.coding:
            info = build_empty_class(self, 'TLseq not in ORF')
            self.utrSeq = 'Sequence unavailable'
            return info

        elif len(TLseq) <= 8:
            info = build_empty_class(self, 'TLseq is too small')
            self.utrSeq = 'Sequence unavailable'
            return info

        elif self.coding not in self.coding or len(self.coding) < len(self.coding):
            info = build_empty_class(self, 'coding not in ORF')
            self.utrSeq = 'Sequence unavailable'
            return info

        else:
            print 'last prob', len(TLseq)






    def MainORFkozak(self):
        five_prime = self.utrSeq[-4:]
        after_atg = self.coding[3:5]
        koz_seq = five_prime + after_atg
        if len(koz_seq) == 6:
            return ScoreKozak(koz_seq)
        else:
            return 0
       


def read_fasta(fp):
        name, seq = None, []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))


classDict = {}
coding_test = []
with open('coding_.txt') as fp:
    for name2, seq2 in read_fasta(fp):
            coding_test.append(name2)
print 'there are', len(coding_test), ' genes in coding.txt file'
cDNAs_, utrs, coding = [], [], []


with open('cDNAs_.txt') as fp11:
    for name2, seq2 in read_fasta(fp11):
        geneID = name2.split('|')[0]
        geneID = re.sub('>', '', geneID)
        transcriptID = name2.split('|')[1]
        gene_name = name2.split('|')[2]
        transcript_type = name2.split('|')[3]
        gene_type = name2.split('|')[4]
        cDNAs_.append(transcriptID)
        geneInfo =  Gene(transcriptID, geneID, gene_name, seq2)
        classDict[geneInfo.transcriptID] = geneInfo
        classDict[transcriptID].transcript_type = transcript_type
        classDict[transcriptID].gene_type = gene_type

print len(classDict), 'classdict after cDNAs', len(cDNAs_)

with open('human_5utr.txt') as fp:
    for name, seq in read_fasta(fp):
        geneID = name2.split('|')[0]
        geneID = re.sub('>', '', geneID)
        transcriptID = name.split('|')[1]
        utrs.append(geneID)
        classDict[transcriptID].utrSeq = seq

      #  classDict[geneInfo.transcriptID] = geneInfo
      #  geneID = Gene(geneID, name)
      #  geneID.utrSeq = seq
      #  classDict[geneID.GeneID] = geneID
  
print len(classDict), 's'
with open('coding_.txt') as fp:
    for name2, seq2 in read_fasta(fp):
            geneID = name2.split('|')[0]
            geneID = re.sub('>', '', geneID)            
            transcriptID = name2.split('|')[1]
            gene_name =  name2.split('|')[2]
            gene_type =  name2.split('|')[4]
            transcript_type =  name2.split('|')[3]
            classDict[transcriptID].coding = seq2
            coding.append(transcriptID)


print 'final dict is', len(classDict)

#for i, t in classDict.iteritems():
  #  print t.coding
#    print t.utrSeq
print len(cDNAs_), len(utrs), len(coding), '[len(cDNAs_), len(utrs), len(coding)]'




print len(classDict)   , 'sss'
################################################################################################
################################################################################################
######################                 set up the classes END                    ###############
################################################################################################
################################################################################################

colNames = [ 'GeneID', 'transcriptID', 'GeneName','TL_Length', 'start', 'stop', 'TL', 'Gene Type', 'transcript_type']


df = pd.DataFrame(columns= colNames)


num = 1
for t, items in classDict.iteritems():


    num = 1+num
    miniDict = {}
    info = items.uORFstart_stop('scrambleEE')
    if info is not None:

        for i in info:

            GeneName = items.GeneName
            transcriptID = items.transcriptID
            transcript_type = items.transcript_type
            gene_type = items.gene_type
            GeneID = items.GeneID
            start = i.start
            stop =  i.stop
            TL = i.TL
            if items.utrSeq == 'Sequence unavailable':
                TL_Length = 'Sequence unavailable'
            else:
                TL_Length = len(items.utrSeq)
            miniDict[num] = GeneID, transcriptID, GeneName,TL_Length, start, stop, TL, gene_type, transcript_type
            dfNew = pd.DataFrame(miniDict, index = colNames).transpose()

            df = df.append(dfNew, ignore_index=True)

    else:

            GeneName = items.GeneName
            gene_type = items.gene_type
            transcriptID = items.transcriptID
            transcript_type = items.transcript_type
            GeneID = items.GeneID
            start= 0
            stop =  0
            TL = 'info is None'
            TL_Length = len(items.utrSeq)
            if items.utrSeq == 'Sequence unavailable':
                TL_Length = 'Sequence unavailable'
            else:
                TL_Length = len(items.utrSeq)
                print 'no uORFs found'
            miniDict[num] = GeneID, transcriptID, GeneName,TL_Length, start, stop, TL, gene_type, transcript_type

            dfNew = pd.DataFrame(miniDict, index = colNames).transpose()
          #  print dfNew
            df = df.append(dfNew, ignore_index=True)

df.to_pickle('human_data3')
df.to_csv('human_data2.csv')

#df = pd.read_pickle('human_data3')
print df
print 'finished..'
        

