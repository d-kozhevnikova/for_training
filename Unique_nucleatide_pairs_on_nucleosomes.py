#!/usr/bin/env python
# coding: utf-8

# In[2]:


import HTSeq
sequences = dict((s.name, s) for s in HTSeq.FastaReader("/home/d_kozhevnikova/genome_s_cer/GCF_000146045.2_R64_genomic.fna"))
sam_reader = HTSeq.SAM_Reader("/home/d_kozhevnikova/sam_files/align49_sorted_by_name.sam")
AA_dict = {a: 0 for a in range(147)}
AT_dict = {a: 0 for a in range(147)}
AC_dict = {a: 0 for a in range(147)}
AG_dict = {a: 0 for a in range(147)}
CA_dict = {a: 0 for a in range(147)}
CC_dict = {a: 0 for a in range(147)}
CG_dict = {a: 0 for a in range(147)}
GA_dict = {a: 0 for a in range(147)}
GC_dict = {a: 0 for a in range(147)}
TA_dict = {a: 0 for a in range(147)}
import re
for alngt in sam_reader.fetch():
    if alngt.aligned:
        if abs(alngt.inferred_insert_size)==147:
            if alngt.mate_aligned == True and (alngt.mate_start.pos < alngt.iv.start)==True:
                #создание <GenomicInterval object 'имя хромосомы, с которой рид', [start,end)>   :
                chrom_name = alngt.iv.chrom
                iv = HTSeq.GenomicInterval( chrom_name, alngt.mate_start.pos, alngt.iv.end, "." )
                str_chrom = sequences[chrom_name].seq
                fragment = str_chrom[alngt.mate_start.pos:alngt.iv.end] #обьект типа <class 'bytes'>
                fragment_seq = str(fragment)[2:149]
                fragment_seq = fragment_seq.upper()
                AA_positions = [m.start() for m in re.finditer('AA', fragment_seq)] 
                #AA_positions - лист с индексами подстрок AA в fragment_seq
                for pos in AA_positions:
                    AA_dict[pos] += 1
                    next_pos = pos + 1
                    AA_dict[next_pos] += 1
                TT_positions = [m.start() for m in re.finditer('TT', fragment_seq)] 
                #TT_positions - лист с индексами подстрок TT в fragment_seq
                for pos in TT_positions:
                    AA_dict[pos] += 1
                    next_pos = pos + 1
                    AA_dict[next_pos] += 1
                AT_positions = [m.start() for m in re.finditer('AT', fragment_seq)] 
                #AT_positions - лист с индексами подстрок AT в fragment_seq
                for pos in AT_positions:
                    AT_dict[pos] += 1
                    next_pos = pos + 1
                    AT_dict[next_pos] += 1
                AC_positions = [m.start() for m in re.finditer('AC', fragment_seq)] 
                #AC_positions - лист с индексами подстрок AC в fragment_seq
                for pos in AC_positions:
                    AC_dict[pos] += 1
                    next_pos = pos + 1
                    AC_dict[next_pos] += 1
                GT_positions = [m.start() for m in re.finditer('GT', fragment_seq)] 
                for pos in GT_positions:
                    AC_dict[pos] += 1
                    next_pos = pos + 1
                    AC_dict[next_pos] += 1
                AG_positions = [m.start() for m in re.finditer('AG', fragment_seq)] 
                for pos in AG_positions:
                    AG_dict[pos] += 1
                    next_pos = pos + 1
                    AG_dict[next_pos] += 1
                CT_positions = [m.start() for m in re.finditer('CT', fragment_seq)] 
                for pos in CT_positions:
                    AG_dict[pos] += 1
                    next_pos = pos + 1
                    AG_dict[next_pos] += 1
                CA_positions = [m.start() for m in re.finditer('CA', fragment_seq)] 
                for pos in CA_positions:
                    CA_dict[pos] += 1
                    next_pos = pos + 1
                    CA_dict[next_pos] += 1
                TG_positions = [m.start() for m in re.finditer('TG', fragment_seq)] 
                for pos in TG_positions:
                    CA_dict[pos] += 1
                    next_pos = pos + 1
                    CA_dict[next_pos] += 1
                CC_positions = [m.start() for m in re.finditer('CC', fragment_seq)] 
                for pos in CC_positions:
                    CC_dict[pos] += 1
                    next_pos = pos + 1
                    CC_dict[next_pos] += 1
                GG_positions = [m.start() for m in re.finditer('GG', fragment_seq)] 
                for pos in GG_positions:
                    CC_dict[pos] += 1
                    next_pos = pos + 1
                    CC_dict[next_pos] += 1
                CG_positions = [m.start() for m in re.finditer('CG', fragment_seq)] 
                for pos in CG_positions:
                    CG_dict[pos] += 1
                    next_pos = pos + 1
                    CG_dict[next_pos] += 1
                GA_positions = [m.start() for m in re.finditer('GA', fragment_seq)] 
                for pos in GA_positions:
                    GA_dict[pos] += 1
                    next_pos = pos + 1
                    GA_dict[next_pos] += 1
                TC_positions = [m.start() for m in re.finditer('TC', fragment_seq)] 
                for pos in TC_positions:
                    GA_dict[pos] += 1
                    next_pos = pos + 1
                    GA_dict[next_pos] += 1
                GC_positions = [m.start() for m in re.finditer('GC', fragment_seq)] 
                for pos in GC_positions:
                    GC_dict[pos] += 1
                    next_pos = pos + 1
                    GC_dict[next_pos] += 1
                TA_positions = [m.start() for m in re.finditer('TA', fragment_seq)] 
                for pos in TA_positions:
                    TA_dict[pos] += 1
                    next_pos = pos + 1
                    TA_dict[next_pos] += 1
                
                
sam_reader2 = HTSeq.SAM_Reader("/home/d_kozhevnikova/sam_files/align50_sorted_by_name.sam")
AT_dict_50 = {a: 0 for a in range(147)}
import re
for alngt in sam_reader2.fetch():
    if alngt.aligned:
        if abs(alngt.inferred_insert_size)==147:
            if alngt.mate_aligned == True and (alngt.mate_start.pos < alngt.iv.start)==True:
                #создание <GenomicInterval object 'имя хромосомы, с которой рид', [start,end)>   :
                chrom_name = alngt.iv.chrom
                iv = HTSeq.GenomicInterval( chrom_name, alngt.mate_start.pos, alngt.iv.end, "." )
                str_chrom = sequences[chrom_name].seq
                fragment = str_chrom[alngt.mate_start.pos:alngt.iv.end] #обьект типа <class 'bytes'>
                fragment_seq = str(fragment)[2:149]
                fragment_seq = fragment_seq.upper()
                AA_positions = [m.start() for m in re.finditer('AA', fragment_seq)] 
                #AA_positions - лист с индексами подстрок AA в fragment_seq
                for pos in AA_positions:
                    AA_dict[pos] += 1
                    next_pos = pos + 1
                    AA_dict[next_pos] += 1
                TT_positions = [m.start() for m in re.finditer('TT', fragment_seq)] 
                #TT_positions - лист с индексами подстрок TT в fragment_seq
                for pos in TT_positions:
                    AA_dict[pos] += 1
                    next_pos = pos + 1
                    AA_dict[next_pos] += 1
                AT_positions = [m.start() for m in re.finditer('AT', fragment_seq)] 
                #AT_positions - лист с индексами подстрок AT в fragment_seq
                for pos in AT_positions:
                    AT_dict[pos] += 1
                    next_pos = pos + 1
                    AT_dict[next_pos] += 1
                AC_positions = [m.start() for m in re.finditer('AC', fragment_seq)] 
                #AC_positions - лист с индексами подстрок AC в fragment_seq
                for pos in AC_positions:
                    AC_dict[pos] += 1
                    next_pos = pos + 1
                    AC_dict[next_pos] += 1
                GT_positions = [m.start() for m in re.finditer('GT', fragment_seq)] 
                for pos in GT_positions:
                    AC_dict[pos] += 1
                    next_pos = pos + 1
                    AC_dict[next_pos] += 1
                AG_positions = [m.start() for m in re.finditer('AG', fragment_seq)] 
                for pos in AG_positions:
                    AG_dict[pos] += 1
                    next_pos = pos + 1
                    AG_dict[next_pos] += 1
                CT_positions = [m.start() for m in re.finditer('CT', fragment_seq)] 
                for pos in CT_positions:
                    AG_dict[pos] += 1
                    next_pos = pos + 1
                    AG_dict[next_pos] += 1
                CA_positions = [m.start() for m in re.finditer('CA', fragment_seq)] 
                for pos in CA_positions:
                    CA_dict[pos] += 1
                    next_pos = pos + 1
                    CA_dict[next_pos] += 1
                TG_positions = [m.start() for m in re.finditer('TG', fragment_seq)] 
                for pos in TG_positions:
                    CA_dict[pos] += 1
                    next_pos = pos + 1
                    CA_dict[next_pos] += 1
                CC_positions = [m.start() for m in re.finditer('CC', fragment_seq)] 
                for pos in CC_positions:
                    CC_dict[pos] += 1
                    next_pos = pos + 1
                    CC_dict[next_pos] += 1
                GG_positions = [m.start() for m in re.finditer('GG', fragment_seq)] 
                for pos in GG_positions:
                    CC_dict[pos] += 1
                    next_pos = pos + 1
                    CC_dict[next_pos] += 1
                CG_positions = [m.start() for m in re.finditer('CG', fragment_seq)] 
                for pos in CG_positions:
                    CG_dict[pos] += 1
                    next_pos = pos + 1
                    CG_dict[next_pos] += 1
                GA_positions = [m.start() for m in re.finditer('GA', fragment_seq)] 
                for pos in GA_positions:
                    GA_dict[pos] += 1
                    next_pos = pos + 1
                    GA_dict[next_pos] += 1
                TC_positions = [m.start() for m in re.finditer('TC', fragment_seq)] 
                for pos in TC_positions:
                    GA_dict[pos] += 1
                    next_pos = pos + 1
                    GA_dict[next_pos] += 1
                GC_positions = [m.start() for m in re.finditer('GC', fragment_seq)] 
                for pos in GC_positions:
                    GC_dict[pos] += 1
                    next_pos = pos + 1
                    GC_dict[next_pos] += 1
                TA_positions = [m.start() for m in re.finditer('TA', fragment_seq)] 
                for pos in TA_positions:
                    TA_dict[pos] += 1
                    next_pos = pos + 1
                    TA_dict[next_pos] += 1
                


# In[3]:


import numpy as np
position = np.array(list(AA_dict.keys())) - 73

valueAA = np.array(list(AA_dict.values()))
valueAT = np.array(list(AT_dict.values()))
valueAC = np.array(list(AC_dict.values()))
valueAG = np.array(list(AG_dict.values()))
valueCA = np.array(list(CA_dict.values()))
valueCC = np.array(list(CC_dict.values()))
valueCG = np.array(list(CG_dict.values()))
valueGA = np.array(list(GA_dict.values()))
valueGC = np.array(list(GC_dict.values()))
valueTA = np.array(list(TA_dict.values()))

val_sumAA = np.sum(valueAA, axis = 0, keepdims = False)
val_sumAT = np.sum(valueAT, axis = 0, keepdims = False)
val_sumAC = np.sum(valueAC, axis = 0, keepdims = False)
val_sumAG = np.sum(valueAG, axis = 0, keepdims = False)
val_sumCA = np.sum(valueCA, axis = 0, keepdims = False)
val_sumCC = np.sum(valueCC, axis = 0, keepdims = False)
val_sumCG = np.sum(valueCG, axis = 0, keepdims = False)
val_sumGA = np.sum(valueGA, axis = 0, keepdims = False)
val_sumGC = np.sum(valueGC, axis = 0, keepdims = False)
val_sumTA = np.sum(valueTA, axis = 0, keepdims = False)

mean_valueAA = val_sumAA/147
mean_valueAT = val_sumAT/147
mean_valueAC = val_sumAC/147
mean_valueAG = val_sumAG/147
mean_valueCA = val_sumCA/147
mean_valueCC = val_sumCC/147
mean_valueCG = val_sumCG/147
mean_valueGA = val_sumGA/147
mean_valueGC = val_sumGC/147
mean_valueTA = val_sumTA/147

value_normalizedAA = valueAA/mean_valueAA
value_normalizedAT = valueAT/mean_valueAA
value_normalizedAC = valueAC/mean_valueAA
value_normalizedAG = valueAG/mean_valueAA
value_normalizedCA = valueCA/mean_valueAA
value_normalizedCC = valueCC/mean_valueAA
value_normalizedCG = valueCG/mean_valueAA
value_normalizedGA = valueGA/mean_valueAA
value_normalizedGC = valueGC/mean_valueAA
value_normalizedTA = valueTA/mean_valueAA

import matplotlib.pyplot as pyplot
pyplot.plot(position, value_normalizedAA)
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of AA dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=10))
pyplot.title(r'AA pairs on nucleosomes: for 3AT dataset', fontsize=13, y=1.05);
pyplot.text(0, 1.0, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [18, 3]
pyplot.show()

pyplot.plot(position, value_normalizedAT)
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of AT dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=10))
pyplot.title(r'AT pairs on nucleosomes: for 3AT dataset', fontsize=13, y=1.05);
pyplot.text(0, 0.6, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [18, 3]
pyplot.show()

pyplot.plot(position, value_normalizedAC)
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of AC dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=10))
pyplot.title(r'AC pairs on nucleosomes: for 3AT dataset', fontsize=13, y=1.05);
pyplot.text(0, 0.9, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [18, 3]
pyplot.show()

pyplot.plot(position, value_normalizedAG)
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of AG dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=10))
pyplot.title(r'AG pairs on nucleosomes: for 3AT dataset', fontsize=13, y=1.05);
pyplot.text(0, 0.8, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [18, 3]
pyplot.show()

pyplot.plot(position, value_normalizedCA)
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of CA dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=10))
pyplot.title(r'CA pairs on nucleosomes: for 3AT dataset', fontsize=13, y=1.05);
pyplot.text(0, 1.0, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [18, 3]
pyplot.show()

pyplot.plot(position, value_normalizedCC)
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of CC dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=10))
pyplot.title(r'CC pairs on nucleosomes: for 3AT dataset', fontsize=13, y=1.05);
pyplot.text(0, 0.53, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [18, 3]
pyplot.show()

pyplot.plot(position, value_normalizedCG)
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of CG dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=10))
pyplot.title(r'CG pairs on nucleosomes: for 3AT dataset', fontsize=13, y=1.05);
pyplot.text(0, 0.225, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [18, 3]
pyplot.show()

pyplot.plot(position, value_normalizedGA)
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of GA dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=10))
pyplot.title(r'GA pairs on nucleosomes: for 3AT dataset', fontsize=13, y=1.05);
pyplot.text(0, 1.0, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [18, 3]
pyplot.show()

pyplot.plot(position, value_normalizedGC)
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of GC dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=10))
pyplot.title(r'GC pairs on nucleosomes: for 3AT dataset', fontsize=13, y=1.05);
pyplot.text(0, 0.3, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [18, 3]
pyplot.show()

pyplot.plot(position, value_normalizedTA)
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of TA dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=10))
pyplot.title(r'TA pairs on nucleosomes: for 3AT dataset', fontsize=13, y=1.05);
pyplot.text(0, 0.55, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [18, 3]
pyplot.show()


# In[20]:


pyplot.plot(position, value_normalizedAT, label = 'AT')
pyplot.plot(position, value_normalizedTA, label = 'TA' )
pyplot.plot(position, value_normalizedAA, label = 'AA' )
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=20))
pyplot.title(r'A/T containing unique nucleatide pairs on nucleosomes', fontsize=13, y=1.05);
pyplot.text(0, 0.1, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [10, 3]
pyplot.legend(loc='upper left')
pyplot.grid(axis = 'x')
pyplot.show()


# In[18]:


o = value_normalizedAT + value_normalizedTA + value_normalizedAA
pyplot.plot(position, o , label = 'AT + TA + AA')
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=20))
pyplot.title(r'A/T containing unique nucleatide pairs on nucleosomes', fontsize=13, y=1.05);
pyplot.text(0, 0.9, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [10, 3]
pyplot.legend(loc='upper left')
pyplot.grid(axis = 'x')
pyplot.show()


# In[25]:


pyplot.plot(position, value_normalizedAT, label = 'AT')
pyplot.plot(position, value_normalizedTA, label = 'TA' )
pyplot.plot(position, value_normalizedAA, label = 'AA' )
pyplot.plot(position, value_normalizedAC, label = 'AС')
pyplot.plot(position, value_normalizedAG, label = 'AG' )
pyplot.plot(position, value_normalizedCA, label = 'CA' )
pyplot.plot(position, value_normalizedCC, label = 'CC')
pyplot.plot(position, value_normalizedCG, label = 'CG' )
pyplot.plot(position, value_normalizedGC, label = 'GC' )
pyplot.plot(position, value_normalizedGA, label = 'GA' )
pyplot.xlabel('Positions on nucleosome (nb)')
pyplot.ylabel('Number of dinucleatides')
pyplot.axis(xmin = -74, xmax = 74, )
pyplot.xticks(np.arange(-70, 75, step=20))
pyplot.title(r'Unique nucleatide pairs on nucleosomes', fontsize=13, y=1.05);
pyplot.text(0, 0.1, 'DYAD',  color = 'red', fontsize=13, horizontalalignment='center')
pyplot.rcParams['figure.figsize'] = [18, 5]
pyplot.legend(loc='upper left')
pyplot.grid(axis = 'x')
pyplot.show()


# In[ ]:




