'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break
def main():
    file="/home/alexanderpetri/Desktop/RAWDATA_PhD1/AccessionPRJEB34849-SIRV/full_length_output.fastq"
    reads = {acc: seq.upper() for acc, (seq, _) in readfq(open(file, 'r'))}
    shortest_seq=1000000000
    for accession,sequence in reads.items():
        print(str(len(sequence)))
        if len(sequence)<10:
            print(sequence)
        if len(sequence) < shortest_seq:
            shortest_seq=len(sequence)
            print("CURRSHORTEST",shortest_seq)
    print("Shortest Sequence:",shortest_seq)

if __name__ == "__main__":
        main()