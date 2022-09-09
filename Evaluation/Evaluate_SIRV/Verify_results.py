import itertools
from consensus import *
from recordclass import recordclass
from IsoformGeneration import align_to_merge
import tempfile
import pickle


def verify():
    max_batch_id=0
    outfolder="out"
    for batchid in range(0,max_batch_id):
        mappingname = "mapping" + str(batchid) + ".txt"
        with open(os.path.join(outfolder, mappingname)) as g:
            for id, reads in itertools.zip_longest(*[g] * 2):
                # print(id, reads)
                inter_id = id.replace('\n', '')
                id = int(inter_id.replace('consensus', ''))
                # print("ID",id)
                reads = reads.replace('\n', '')
                reads = reads.replace('[', '')
                reads = reads.replace(']', '')
                # reads=reads.replace('\'','')
                reads = reads.replace("'", "")
                readslist = reads.split(", ")
                # print(readslist)
                # print(len(readslist))
                batch_mappings[id] = readslist