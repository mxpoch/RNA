# imports
import numpy as np
import pandas as pd
import math 
import pickle
import re
import json
import sys
import os
sys.setrecursionlimit(10000)

##Database Formatters
# Accession | Name | Host | Additional Details | Nucleotide completeness |
def rna_metadata_extractor(desc: str) -> list:
    # splits input string by |
    metadata = desc.split('|')
    # moves 'Additional Details' to 3rd index
    temp = metadata.pop(1)
    metadata.insert(3, temp)
    return metadata

def import_rna(rna_fasta_path: str) -> pd.DataFrame():
    cols = ['accession','name','host','additional_details','nucleotide_completeness','length','seq']
    # created dataframe of arbitrary size 
    nucDB = pd.DataFrame(index=range(1_000_000), columns=cols)
    with open(corona_nucleotides) as corona_import:
        # creates generator object for each entry in .fasta file (~2.0GB)
        crowns = SeqIO.parse(corona_import, 'fasta')
        i = 0
        for virus in crowns:
            # extracts metadata
            row_buffer = rna_metadata_extractor(virus.description)
            # creates sequence length column, useful for filtering 
            length = len(virus.seq)
            row_buffer.append(length)
            # appends actual RNA sequence
            row_buffer.append(str(virus.seq))
            nucDB.loc[i] = row_buffer
            # print(i)
            i += 1
    # (deallocates?) empty rows
    nucDB.dropna(how='all', axis=0, inplace=True)
    return nucDB
# time taken: 27m 46s


## Global Functions
# converts subsquence ranges to sequences 
def convert_to_sequence(origin_seq:str, ranges:list) -> list:
        seqs = []
        for x in ranges:
            seqs.append([x[0], origin_seq[x[0]:x[1]]])
        return seqs
# converts one subsequence range to sequence 
def convert_one_to_sequence(seq:str, range:list) -> str: return seq[range[0]:range[1]]
# verifier
# the most important function here.
def verifier(vs1:list, vs2:list) -> list:
    nlist = []
    for x in vs1:
        flag = False
        for y in vs2:
            if re.search(x[1], y[1]):
                flag = True
                break
        if flag == False:
            nlist.append(x)
    # {key} : [locations in seq1 (validation strand)]
    no_dupli = {}
    for seg in nlist:
        # print(ss_dict.keys())
        if seg[1] not in no_dupli.keys(): no_dupli[seg[1]] = [seg[0]]
        else: no_dupli[seg[1]].append(seg[0])
    # prints lengths of sequence being compared against, the number of unique missing sequences
    # and the missing sequences themselves
    print(len(vs1), len(no_dupli.keys()), no_dupli)

## Algorithm Classes
class SuffixTreePreProcess:
    def __init__(self, TARGET:pd.DataFrame, alphabet:list):
        self.DB = TARGET
        self.codons = self.create_codons(alphabet)

    def create_codons_worker(self, clist: list, length:int, alphabet: list, ret:list):
        for a in alphabet:
            if length == 3:
                ret.append(''.join(clist))
                return
            clist.append(a)
            self.create_codons_worker(clist, length+1, alphabet, ret)
            clist.pop(-1)
    
    def create_codons(self, alphabet:list) -> None:
        temp = []
        self.create_codons_worker([], 0, alphabet, temp)
        codons = [[-1, x] for x in temp]
        return codons
    
    # generates sources based on the list of indices (selected:list) of the pandas dataframe 
    def generate_sources(self, selected:list) -> tuple:
        TSEQ = [self.DB.iloc[x] for x in selected]
        SEQ = [X['seq'] for X in TSEQ]
        return SEQ, TSEQ
    
    def fill_indices(self, seq:str, cddict: dict):
        for i in range(len(seq)-2):
            cframe = seq[i:i+3]
            if cframe not in cddict.keys(): continue
            cddict[cframe].append(i)
        return cddict
    
    def create_dictionary(self, keys: list) -> dict:
        frame_dict = {}
        for k in keys:
            frame_dict[k[1]] = []
        return frame_dict

    def generate_returns(self, SEQ) -> tuple:
        RET = [[] for x in range(len(SEQ))] # output array
        SS = [[] for x in range(len(SEQ))] # temp working array
        CT = [self.fill_indices(x, self.create_dictionary(self.codons)) for x in SEQ] # codon tables
        
        return RET, SS, CT

class SuffixTreeSearch:
    def __init__(self, SEQ, RET, SS, CT):
        self.SEQ = SEQ
        self.RET = RET
        self.SS = SS
        self.CT = CT
        self.stack = []

    def binary_search(self, search:int, li:list) -> int:
        # single value edge-case
        if len(li) == 0:
            return -1
        if len(li) == 1:
            if li[0] == search: return 0
            else: return -1

        top = len(li)-1
        bottom = 0
        i = math.ceil((top-bottom)/2)
        icounter = 0
        prev = 0 

        # print(len(li), i)
        while search != li[i]:
            # small region edge-case
            if top-bottom <= 3:
                for i in range(bottom, top):
                    if li[i] == search:
                        return i
                return -1
            
            if search > li[i]:
                bottom = i
                i = math.ceil((top-i)/2) + bottom
            elif search < li[i]:
                top = i
                i = math.ceil((i-bottom)/2) + bottom
            if i == prev:
                icounter += 1
                if icounter > 1:
                    return -1
            prev = i 
        return i

    def check_base(self) -> bool:
        for ss in self.SS:
            for x in ss:
                if x[2] == 0:
                    return True
        return False

    def check_empty(self) -> bool:
        for ss in self.SS:
            if len(ss) > 0:
                return False
        return True

    def consolidate_keys(self, recur_depth:int) -> list:
        TEMP_KEYS = [[] for x in range(len(self.SS))]

        i = 0 # sync 
        for ss in self.SS:
            for ss_e in ss:
                # finds appropriate length sequences
                if ss_e[2] == recur_depth and ss_e[1]+6 <= len(self.SEQ[i]):
                        # format: (hash of original sequence, next codon)
                        hkey = (ss_e[3], self.SEQ[i][ss_e[1]+3:ss_e[1]+6])
                        TEMP_KEYS[i].append(hkey)
            i += 1
        return TEMP_KEYS

    def check_keys(self, KEYS:list) -> list:
        validk = []
        # takes the first set of keys, and compares against the entire set
        # only keys common across all are returned
        for validator in KEYS[0]:
            for kcheck in KEYS[1:]:
                if validator in kcheck and validator not in validk:
                    validk.append(validator)
                    break
        return validk

    def valid(self, key:str) -> bool:
        for dic in self.CT:
            if len(dic[key]) > 0:
                return True
        return False

    def base_key(self, ss:list, key:list, ct:dict):
        for indx in ct[key]:
            # ss format = index start, index end, recur_depth/length, hashkey
            ss.append([indx, indx, 0, hash(key)])
    
    def shadow_key_subroutine(self, SS:list, CT:list, SEQ:list):
        # print('shadowkey')
        # 1. Find reverse_key
        # inits reverse key array
        KEYS = [[] for x in range(len(SS))]
        i = 0 # to sync SS, CT, and SEQ
        for ss in SS:
            for ss_e in ss:  
                if ss_e[0]-3 > 0:
                    current = SEQ[i][ss_e[0]-3:ss_e[0]+3]
                    KEYS[i].append(current)  
            i += 1

        # 2. Validate reverse_key
        valid_keys = self.check_keys(KEYS)

        # 3. Grow reverse sequence if valid 
        # DOES NOT GROW recur_depth (ss[2])
        j = 0 # sync to SS, CT, and SEQ
        for ss in SS:
            # subsequence lists
            for ss_e in ss:
                # actual subsequences
                for valid in valid_keys:
                    # scanning through all 
                    if valid == SEQ[j][ss_e[0]-3:ss_e[0]+3] and ss_e[0]-3 > 0:
                        ss_e[0] -= 3
                        ss_e[3] = hash(SEQ[j][ss_e[0] : ss_e[0]+6])
                        break
            j += 1

    def refresh_ss(self, SS:list, RET:list, recur_depth:int):
        for r in range(len(SS)):
            # print(SS[r])
            i = 0 
            while i < len(SS[r]):
                # format: SS[ss][ss_e][indx]
                if SS[r][i][2] == recur_depth:
                    SS[r][i][1] = SS[r][i][1] + 3 # corrects range
                    RET[r].append(SS[r][i][:2])
                    del SS[r][i]
                else:
                    i += 1 
    
    def grow_seq(self, ss:list, key:list, seq:str, recur_depth:int):
        for i in range(len(ss)):
            if ss[i][3] == key[0] and seq[ss[i][1]+3:ss[i][1]+6] == key[1] and ss[i][2] == recur_depth-1:
                ss[i][1] = ss[i][1] + 3
                ss[i][2] = ss[i][2] + 1
                ss[i][3] = hash(seq[ss[i][0]:ss[i][1]+3])
            else: continue
    
    def convertt(self, ss:list) -> list:
        c = []
        for s in ss:
            current = s.copy()
            current[1] += 3
            c.append(current)
        return c

    def count_duplicates(self, raw_comm:list) -> dict:
        ret = {}
        for x in raw_comm:
            if tuple(x) in ret.keys(): continue
            ccount = raw_comm.count(x)
            if ccount > 1: ret[tuple(x)] = ccount
        return ret

    def find_common_worker(self, keys:list, recur_depth:int): 
        # [initial keys, recursion depth, and loop counter (must be < len(keys))]
        # self.stack[len(self.stack)-1]... => last entry in the FIFO stack
        self.stack.append([keys, recur_depth, 0])

        while len(self.stack) > 0: 
            # print(self.stack[len(self.stack)-1][0], self.stack[len(self.stack)-1][2], self.stack[len(self.stack)-1])
            # key = [origin sequence hash, adjusted hash, key]
            key = self.stack[len(self.stack)-1][0][self.stack[len(self.stack)-1][2]]
            # incrementing the loop counter
            if len(self.stack[len(self.stack)-1][0]) > self.stack[len(self.stack)-1][2]:
                self.stack[len(self.stack)-1][2] += 1
            # print('Current keyset:', [k for k in keys])
            if self.stack[len(self.stack)-1][1] == 0:
                print('rootkey', key)
            # clears out keys for next round
            # if SS contains bases and the recursion depth is 0
            if self.stack[len(self.stack)-1][1] == 0: 
                for ss in self.SS: ss.clear()
            # ensures the key arent empty in both codon tables (rare chance, but it may happen)
            if self.valid(key[1]):
                # if all sets are empty and recur_depth is at 0
                # create base keys from the current key. 
                if self.stack[len(self.stack)-1][1] == 0:
                    i = 0
                    for ss in self.SS: 
                        # base keys are (start, stop-3, recur_depth, hash)
                        self.base_key(ss, key[1], self.CT[i]) 
                        i += 1
                    self.shadow_key_subroutine(self.SS, self.CT, self.SEQ)

                # grows sequences 
                else:
                    i = 0 
                    for ss in self.SS: 
                        self.grow_seq(ss, key, self.SEQ[i], self.stack[len(self.stack)-1][1])
                        i += 1 

                # finds new keys
                TEMP_KEYS = self.consolidate_keys(self.stack[len(self.stack)-1][1])
                # filters out new keys
                search_keys = self.check_keys(TEMP_KEYS)
                
                # iterative management
                # appending new iterator
                if len(self.stack[len(self.stack)-1][0]) > 0 and len(search_keys) > 0:
                    self.stack.append([search_keys, self.stack[len(self.stack)-1][1]+1, 0])
            # makes sure there is at least 1 common subsequence in this depth of the recursion tree
            self.refresh_ss(self.SS, self.RET, self.stack[len(self.stack)-1][1])
            # removing expired iterators
            # while stack exists and loop counter is greater than the number of keys
            while len(self.stack) > 0 and len(self.stack[len(self.stack)-1][0]) <= self.stack[len(self.stack)-1][2]:
                del self.stack[len(self.stack)-1]       

    def find_common(self, keys:list):
        self.find_common_worker(keys, 0)

    def cleanse_duplicates(self):
        j = 0
        for ret in self.RET:
            tstring = convert_to_sequence(self.SEQ[j], ret)
            i = 0
            while i < len(tstring):
                ss = tstring[i]

                while tstring.count(ss) > 1:
                    cindx = tstring.index(ss)
                    del tstring[cindx]
                    del ret[cindx]
                
                i += 1
            j += 1

class SuffixTreePostProcess:
    def __init__(self, RET:list, SEQ:list, TSEQ:list):
        self.IMPORT = RET
        self.SEQ = SEQ
        self.TSEQ = TSEQ

    def convert_to_ss(self, converted:list, index_dict:dict):
        index_in_seq = 0
        for subseq in converted:
            if subseq[1] not in index_dict.keys(): # replace for hash keys 
                index_dict[subseq[1]] = {'subseq': [(subseq[0], index_in_seq)], 'merged': False}
            else:
                index_dict[subseq[1]]['subseq'].append((subseq[0], index_in_seq))
            index_in_seq += 1

    def find_like(self, left:tuple, right:tuple, index_dict:dict, overlap_amount:int) -> list:
        # checking for edge cases
        if left[1] not in index_dict.keys() or right[1] not in index_dict.keys():
            # print('a', end=' ')
            return [(-1, -1, -1, -1)]
        if overlap_amount < 0:
            # print('b', end=' ')
            return [(-1, -1, -1, -1)]

        likewise = []
        for like_left in index_dict[left[1]]['subseq']:
            for like_right in index_dict[right[1]]['subseq']:
                compliment_overlap = like_left[0]+len(left[1]) - like_right[0]
                
                if compliment_overlap == overlap_amount:
                    # stop determines the longer subsequences (end index)
                    stop = max(like_right[0]+len(right[1]), like_left[0]+len(left[1]))
                    #                start         stop                     index in seq -left  index in seq -right
                    likewise.append((like_left[0], stop, like_left[1], like_right[1]))
                    break
        if len(likewise) == 0:
            return [(-1,-1,-1,-1)]
        return likewise

    def valid_merge_indices(self, MERGE_INDICES:list) -> bool:
        # entire merge_indices array
        for MI in MERGE_INDICES:
            # entries in entire array
            for mi in MI:
                if mi[0] == -1:
                    return False
        return True

    def update_offsets(self, right_index:int, index_offsets:dict):
        if right_index in index_offsets.keys(): index_offsets[right_index] += 1
        else: index_offsets[right_index] = 1

    def cleanse_index_dict(self, left:list, right:list, m_and_a:tuple, index_dict:dict):
        index_dict[left[1]]['subseq'] = [x for x in index_dict[left[1]]['subseq'] if x[1] != m_and_a[2]]
        index_dict[right[1]]['subseq'] = [x for x in index_dict[right[1]]['subseq'] if x[1] != m_and_a[3]]
        if len(index_dict[left[1]]['subseq']) == 0: del index_dict[left[1]]
        if len(index_dict[right[1]]['subseq']) == 0: del index_dict[right[1]]

    def remove_ss_inplace(self, left:list, right:list, index_dict:list, index_offsets:list, merge_indices:list, raw_comm:list, seq:str):
        # have to implement index offsets for lowest ranges...
        for m_and_a in merge_indices:
            # m_and_a = (start, stop, index position of left, index position of right)
            # deletes index left
            
            index_offset_left = 0 
            #                              just >, not >=
            select_indices_left = list(filter(lambda x: x < m_and_a[2], index_offsets.keys()))
            for s in select_indices_left: index_offset_left += index_offsets[s]
            # deletes merged left ss from raw_comm 
            del raw_comm[m_and_a[2]-index_offset_left]
            
            # deletes right index
            index_offset_right = 0
            select_indices_right = list(filter(lambda x: x < m_and_a[3], index_offsets.keys()))
            for s in select_indices_right: index_offset_right += index_offsets[s]   
            del raw_comm[m_and_a[3]-index_offset_right-1]

            raw_comm.insert(m_and_a[3]-index_offset_right-1, [m_and_a[0], m_and_a[1]])
            
            self.update_offsets(m_and_a[3], index_offsets)
            # removes old index positions from index_dict so that they are not considered in future merge events
            self.cleanse_index_dict(left, right, m_and_a, index_dict)

            # updates index_dict
            # inserts new range into index_dict
            current_key = seq[ m_and_a[0]:m_and_a[1] ]
            if current_key not in index_dict.keys(): 
                index_dict[current_key] = {'subseq':[(m_and_a[0], m_and_a[2])], 'merged':True}
            else: 
                index_dict[current_key]['subseq'].append((m_and_a[0], m_and_a[2]))
                # index_dict[current_key]['merged'] = True # may be redundant -- test

    def merge_raw(self) -> tuple:
        # initalizing strands 
        RAW_COMM = [sorted(x) for x in self.IMPORT]
        # print(RAW_COMM[0], len(RAW_COMM[0]))
        # initializing index dictionaries
        CONVERTED = [convert_to_sequence(self.SEQ[i], RAW_COMM[i]) for i in range(len(RAW_COMM))]
        INDEX_DICT = [{} for x in CONVERTED]

        # converts ss to index dict for easy traversal
        # formatting = {'Sequence': (index in seqeunce, original index in raw_comm list)}
        for i in range(len(CONVERTED)): self.convert_to_ss(CONVERTED[i], INDEX_DICT[i])
        # print(INDEX_DICT[0])

        i = 1
        INDEX_OFFSETS = [{} for x in range(len(INDEX_DICT))]  # (index) --> always offsets by 1
        while i < len(RAW_COMM[0]): # follows only validation strand 
            # traverses validation strand in pairs
            # creates left and right pair for merging
            left = [RAW_COMM[0][i-1][0], convert_one_to_sequence(self.SEQ[0], RAW_COMM[0][i-1])]
            right = [RAW_COMM[0][i][0], convert_one_to_sequence(self.SEQ[0], RAW_COMM[0][i])]

            # calculates overlap amount 
            overlap_amount = RAW_COMM[0][i-1][1] - RAW_COMM[0][i][0]
            # print(RAW_COMM[0][i-1][1], RAW_COMM[0][i][0])
            if overlap_amount < 0:
                i += 1
                continue
            # print(overlap_amount)

            # finds all other identically overlapping sequences 
            MERGE_INDICES = [self.find_like(left, right, index_dict, overlap_amount) for index_dict in INDEX_DICT]
            # print(MERGE_INDICES, CONVERTED[1][982], CONVERTED[0][982])
            # should return list of 4-tuples containing (start, stop, index position of left, index position of right)
            
            # validator for checking with the other index dictionaries
            if not self.valid_merge_indices(MERGE_INDICES): 
                i += 1
                continue
            # uses index_offset to delete seqeunces being merged
            # inserts new, merged range 
            # have to implement index offsets for lowest ranges...
            for i in range(len(RAW_COMM)): self.remove_ss_inplace(left, right, INDEX_DICT[i], INDEX_OFFSETS[i], MERGE_INDICES[i], RAW_COMM[i], self.SEQ[i])   

        for seq in INDEX_DICT:
            for ky in seq.keys():
                seq[ky]['subseq'] = [x[0] for x in seq[ky]['subseq']]
                    
        return RAW_COMM, INDEX_DICT

    def remove_inexact(self, INDEX_DICT:list) -> dict:
        nlist = []
        for vk in INDEX_DICT[0].keys():
            valid = True
            for ids in INDEX_DICT[1:]:
                if vk not in ids.keys(): 
                    valid = False
                    break
            if valid or INDEX_DICT[0][vk]['merged']: nlist.append(vk)
        
        RET = [{} for sequence in INDEX_DICT]
        i = 0 
        for seq in INDEX_DICT:
            for valid in nlist:
                RET[i][valid] = seq[valid]
            i += 1
        
        return RET

    def ss_to_json(self, EXACT:list, filename:str) -> dict:
        jsonDICT = {}
        # creating metadata
        meta = {}
        i = 0
        for exact in EXACT:
            current = {}
            current['len'] = len(self.TSEQ[i]['seq'])
            current['hash'] = self.TSEQ[i]['accession']
            meta[i+1] = current
            i += 1
        jsonDICT['sequences'] = meta
        # creating actual data
        j = 0 
        for ss in EXACT[0].keys():
            current = {}
            current['ss_len'] = len(ss)
            current['ss_hash'] = hash(ss)
            k = 0
            for exact in EXACT: 
                name = self.TSEQ[k]['accession']
                current[name] = exact[ss]
                k += 1
            jsonDICT[f'ss{j}'] = current
            j += 1

        with open(f'/content/drive/MyDrive/{filename}.json', 'w') as viz:
            json.dump(jsonDICT, viz, indent=4)