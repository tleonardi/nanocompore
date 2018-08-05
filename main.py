import logging
from time import sleep
import queue 
import threading
import mmap
import pysam
from collections import defaultdict


logging.basicConfig(level=logging.DEBUG, format='%(threadName)s %(message)s')
logger = logging.getLogger('simple_example')
#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

class events_combinator(object):
    """ Produce usefule results """
    
    def __init__(self, whitelist):
        self.__main_queue = queue.Queue(maxsize=20)
        self.results=[] 
        self.data_f1=dict()
        self.data_f2=dict()
        self.processed=set()
        self.n_consumer_threads=4
        self.lock = threading.Lock()
        self.tx_whitelist=whitelist
        self.file1="/mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_wt/nanopolish/events.txt"
        #self.file1="events.txt"
        #self.file2="events2.txt"
        self.file2="/mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_kd/nanopolish/events.txt"
        
        self.reader1 = threading.Thread(name="reader1", target=self.__parse_events_file, args=(self.file1, self.data_f1,))
        self.reader2 = threading.Thread(name="reader2", target=self.__parse_events_file, args=(self.file2, self.data_f2,))
        
        consumers=[]
        for i in range(self.n_consumer_threads):
            consumers.append(threading.Thread(name="consumer%s"%i, target=self.__consumer))
        
        self.reader1.start()
        self.reader2.start()
        
        for c in consumers:
            c.start()
        
        self.reader1.join()
        self.reader2.join()
        self.__main_queue.join()
        for c in consumers:
            c.join()
        for r in self.results:
            print(r)


    @staticmethod
    def mmap_readline(file, sep="\t"):
        """ Iterator for opening files as mem maps """
        with open(file, "r+b") as f:
            m = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
            while True:
                line=m.readline()
                if line == b'': break
                yield line.decode().split(sep)[0:7]
     
    def __parse_events_file(self, file, data):
        """ Parse an events file line by line and extract lines
        corresponding to the same transcripts. When a transcript
        is read entirely, add its name to the __main_queue and the
        the data to data. """
        tx=""
        res=[]
        for line in self.mmap_readline(file):
            if(line[0] in self.tx_whitelist):
                if line[0]==tx:
                    res.append(line)
                    tx=line[0]
                else:
                    if len(res)!=0:
                        # Might be worth making res an object of a dedicated class
                        data[tx]=res
                        self.__main_queue.put(tx)
                    logger.debug("Finished processing (%s)" % (tx))
                    res=[line]
                    tx=line[0]
        if len(res)!=0:
            data[tx]=res
            self.__main_queue.put(tx)
            logger.debug("Finished processing (%s)" % (tx))

    def __consumer(self):
        while True:
            if self.__main_queue.empty() and not self.reader1.is_alive() and not self.reader2.is_alive():
                logger.debug("Queue is empty and readers have finished, terminating")
                break
            # This timeout prevents threads to get stuck in case the 
            # queue gets empty between the check above and the get
            try:
                tx = self.__main_queue.get(timeout=3)
            except queue.Empty:
                continue

            # If the data is not ready 
            if tx not in self.data_f1.keys() or tx not in self.data_f2.keys():
                # If one of the readers is still alive return tx to the queue
                if self.reader1.is_alive() or self.reader2.is_alive():
                    self.__main_queue.put(tx)
                    logger.debug("Returning (%s) to Q. Q has (%s) items" % (tx, self.__main_queue.qsize()))
                # Else, if either of the readers has finished, remove the task from queue
                else:
                    logger.debug("Discarding (%s) because readers have terminated and data not present" % (tx))
                self.__main_queue.task_done()

            # If the data is ready   
            else: 
                logger.debug("Picked (%s) from Q. Q has (%s) items." % (tx, self.__main_queue.qsize()))
                # Acquire lock to make sure nothing is written to processed after we check
                self.lock.acquire()
                if tx not in self.processed:
                    self.processed.add(tx)
                    self.lock.release()
                    self.results.append(self.process_tx(self.data_f1[tx], self.data_f2[tx]))
                else:
                    self.lock.release()
                self.__main_queue.task_done()

    @staticmethod
    def process_tx(data1, data2):
        if data1[0][0] != data2[0][0]:
            print("Error")
        else: tx=data1[0][0]
        logger.info("Processed %s" % (tx))
        return("Tx (%s) done. Len f1: (%s) - Len f2: (%s)" % (tx, len(data1), len(data2)))
    
class alignements(object):
    def __init__(self, bam):
        self.__samfile = pysam.AlignmentFile(bam, "rb")

    def tx_coverage_filter(self, n=10):
        """ Return the name of transcripts with more than n primary alignemnts """
        tx_sam_counts=defaultdict(int)
        for record in self.__samfile:
            if not record.flag & (256+4):
                tx_sam_counts[record.reference_name]+= 1   
        return set( [ k for (k,v) in tx_sam_counts.items() if v>n ])



if __name__ == '__main__':
    bam1="/mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_wt/minimap/aln_transcriptome.sorted.bam"
    bam2="/mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_kd/minimap/aln_transcriptome.sorted.bam"
    aln_f1=alignements(bam1).tx_coverage_filter(n=10)
    aln_f2=alignements(bam2).tx_coverage_filter(n=10)
    whitelist=set.intersection(aln_f1, aln_f2)
    d = events_combinator(whitelist=whitelist)
