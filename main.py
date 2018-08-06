import logging
import queue 
import threading
import mmap
import pysam
import argparse
import pickle
import os
from tqdm import tqdm
from pathlib import Path
from collections import defaultdict
from os import makedirs


logging.basicConfig(level=logging.DEBUG, format='%(threadName)s %(message)s')
logger = logging.getLogger('nanocompore')
logLevel=logging.WARNING

class nanocompore(object):
    """ Produce usefule results """
    
    def __init__(self, file1, file2, outfolder=None, nthreads=4, whitelist=None):
        """ Main routine that starts readers and consumers 
            file1: path to file1
            file2: path to file2
            outfolder: path to folder where to save output
            nthreads: number of consumer threads
            whitelist: list of transcripts to process. If None (default) process all transcripts
        """

        # Check that input files exist
        for f in (file1, file2):
            if not Path(f).exists():
                raise nanocomporeError("%s not found)" % f)
        self.file1=file1
        self.file2=file2


        # Check that output folder doesn't exist and create it
        if Path(outfolder).exists():
            raise nanocomporeError("The output folder specified already exists")
        else:
            try: 
                makedirs(outfolder)
            except:
                raise nanocomporeError("Error creating output folder %s" % outfolder)
            self.__outfolder=outfolder

        # Check thread number is valid
        try:
            self.n_consumer_threads=int(nthreads)
        except:
            raise nanocomporeError("Number of threads not valid")

        # Check whitelist is valid
        if not isinstance(whitelist, set) and whitelist is not None:
            raise nanocomporeError("Whitelist not valid")
        else:
            self.tx_whitelist=whitelist


    def process(self, checkpoint=True):
        max_qsize=5000
        # Main processing queue
        self.__main_queue = queue.Queue(maxsize=max_qsize)
        # Threading lock
        self.__lock = threading.Lock()
        # Final results
        self.results=[] 
        # Data read by the readers
        self.data_f1=dict()
        self.data_f2=dict()
        # Set containg names of already processed transcripts
        self.processed=set()


        self.__reader1 = threading.Thread(name="reader1", target=self.__parse_events_file, args=(self.file1, self.data_f1, 1))
        self.__reader2 = threading.Thread(name="reader2", target=self.__parse_events_file, args=(self.file2, self.data_f2, 2))

        consumers=[]
        for i in range(self.n_consumer_threads):
            consumers.append(threading.Thread(name="consumer%s"%i, target=self.__consumer))
        
     
        self.__reader1.start()
        self.__reader2.start()
        
        for c in consumers:
            c.start()

        # Wait for the readers to complete    
        self.__reader1.join()
        self.__reader2.join()
        # Wait for the queue to be empty
        logger.warning("Wainting for processing queue to be empty")
        self.__main_queue.join()
        # Wait for the consumers to complete
        for c in consumers:
            c.join()
        # Print results
        if checkpoint:
            with open(self.__outfolder+'/results.p', 'wb') as pfile:
                pickle.dump(self.results, pfile)

    @staticmethod
    def mmap_readline(file, sep="\t"):
        """ Iterator for opening files as mem maps """
        with open(file, "r+b") as f:
            m = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
            while True:
                line=m.readline()
                if line == b'': break
                yield line.decode().split(sep)
     
    def __parse_events_file(self, file, data, barPos=1):
        """ Parse an events file line by line and extract lines
        corresponding to the same transcripts. When a transcript
        is read entirely, add its name to the __main_queue and the
        the data to data. """
        tx=""
        block=[]
        # Reading the whole file to count lines is too expensive. 148.52bytes is the average line size in the events files
        bar=tqdm(total=int(os.stat(file).st_size/148.52), desc="File%s progress" % barPos, position=barPos, unit_scale=True, disable=(logLevel in [logging.DEBUG, logging.INFO]))
        for line in self.mmap_readline(file):
            bar.update()
            if(line[0] != "contig" and line[0] in self.tx_whitelist):
                if line[0]==tx:
                    block.append(line)
                    tx=line[0]
                else:
                    if len(block)!=0:
                        data[tx]=block
                        self.__main_queue.put(tx)
                    logger.debug("Finished processing (%s)" % (tx))
                    block=[line]
                    tx=line[0]
        if len(block)!=0:
            data[tx]=block
            self.__main_queue.put(tx)
            logger.debug("Finished processing (%s)" % (tx))
        bar.close()

    def __consumer(self):
        while True:
            if self.__main_queue.empty() and not self.__reader1.is_alive() and not self.__reader2.is_alive():
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
                if self.__reader1.is_alive() or self.__reader2.is_alive():
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
                self.__lock.acquire()
                if tx not in self.processed:
                    self.processed.add(tx)
                    self.__lock.release()
                    self.results.append(self.process_tx(self.data_f1[tx], self.data_f2[tx]))
                else:
                    self.__lock.release()
                self.__main_queue.task_done()

    @staticmethod
    def process_tx(data1, data2):
        if data1[0][0] != data2[0][0]:
            print("Error")
        else: tx=data1[0][0]
        logger.info("Processed %s" % (tx))
        return([tx, [data1, data2]])
    
class alignements(object):
    def __init__(self, bam):
        # Check that bam file exists
        if not Path(bam).exists():
            raise nanocomporeError("Bam file %s not found)" % bam)
        self.__samfile = pysam.AlignmentFile(bam, "rb")

    def tx_coverage_filter(self, n=10):
        """ Return the name of transcripts with more than n primary alignemnts """
        try:
            n = int(n)
        except:
            raise nanocomporeError("Read coverage argument for whitelisting is not valid")
        tx_sam_counts=defaultdict(int)
        for record in self.__samfile:
            # Filter alignmenets without "Secondary alignment" or "Not mapped" flags
            if not record.flag & (256+4):
                tx_sam_counts[record.reference_name]+= 1   
        return set( [ k for (k,v) in tx_sam_counts.items() if v>n ])



def main():
    parser = argparse.ArgumentParser(description="""Find differences in two nanopolish events files""")
    parser.add_argument('--version', '-v', action='version', version="0.1a")
    parser.add_argument("--file1", required=True,
        help="Path to the Nanopolish events file for sample 1")
    parser.add_argument("--file2", required=True,
        help="Path to the Nanopolish events file for sample 2")
    parser.add_argument("-o", "--outfolder", required=True,
        help="Path to a directory where to store the results. Must not exist")
    parser.add_argument("-n", type=int, default=6, required=False,
        help="Number of threads (two are used for reading files, all the other for processing in parallel).")
    parser.add_argument("--bam1", required=True,
        help="Path to the Minimap2 BAM file for sample 1")
    parser.add_argument("--bam2", required=True,
        help="Path to the Minimap2 BAM file for sample 2")
    parser.add_argument("--mincoverage", type=int, default=0, required=False,
        help="Minimum number of reads covering a transcript (in each sample) for it to be considered. Set to 0 to consider all transcripts.")
    parser.add_argument("--logLevel", default="warning",
        help="Set the log level. Valida values: warning, info, debug.")
    args = parser.parse_args()
    logLevel=args.logLevel
    if args.logLevel == "debug":
        logLevel=logging.DEBUG
    elif args.logLevel == "warning":
        logLevel=logging.WARNING
    elif args.logLevel == "info":
        logLevel=logging.INFO
    else:
        raise nanocomporeError("Log level not valid")
    
    logger.setLevel(logLevel)

    logger.warning("Reading BAM alignments")
    aln_f1=alignements(args.bam1).tx_coverage_filter(n=args.mincoverage)
    aln_f2=alignements(args.bam2).tx_coverage_filter(n=args.mincoverage)
    whitelist=set.intersection(aln_f1, aln_f2)
    logger.warning("Keeping %s transcripts for processing" % len(whitelist))
    n = nanocompore(file1=args.file1, file2=args.file2, nthreads=args.n, whitelist=whitelist, outfolder=args.outfolder)
    logger.warning("Starting multi-threaded processing of events files")
    n.process()
    logger.warning("All completed")
    
class nanocomporeError(Exception):
    """ Basic exception class for nanocompore module """
    pass

if __name__ == '__main__':
    main()
