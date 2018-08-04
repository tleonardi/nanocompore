import logging
from time import sleep
from queue import Queue
import threading
import mmap


logging.basicConfig(level=logging.DEBUG, format='%(threadName)s %(message)s')
logger = logging.getLogger('simple_example')
logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)


class events_combinator(object):
    """ Produce usefule results """
    
    def __init__(self):
        self.__main_queue = Queue(maxsize=20)
        self.results=[] 
        self.data_f1=dict()
        self.data_f2=dict()
        self.processed=set()
        self.n_consumer_threads=3
        
        #self.file1="/mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_wt/nanopolish/events.txt"
        self.file1="events.txt"
        self.file2="events2.txt"
        #self.file2="/mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_wt/nanopolish/events.txt"
        
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
        for r in self.results:
            print(r)
        for _ in consumers:
            self.__main_queue.put(None)


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
            if line[0]==tx:
                res.append(line)
                tx=line[0]
                logger.debug("Processing (%s)" % (tx))
            else:
                if len(res)!=0:
                    # Might be worth making res an object of a dedicated class
                    data[tx]=res
                    self.__main_queue.put(tx)
                logger.debug("Finished processing (%s)" % (tx))
                res=[line]
                tx=line[0]
        data[tx]=res
        self.__main_queue.put(tx)
        logger.debug("Finished processing (%s)" % (tx))

    def __consumer(self):
        while True:
            # Fetch a tx from the queue
            tx =  self.__main_queue.get()
            logger.debug("Picked (%s) from Q. Q has (%s) items." % (tx, self.__main_queue.qsize()))
            # Terminate consumer is we reach
            # bottom of the queue
            if tx is None:
                break
            # Check data is in data_1 and data_2, if not put back in __main_queue
            if tx not in self.processed:
                if tx in self.data_f1.keys() and tx in self.data_f2.keys():
                    self.results.append(self.process_tx(self.data_f1[tx], self.data_f2[tx]))
                    self.processed.add(tx)
                elif self.reader1.is_alive() or self.reader2.is_alive():
                    self.__main_queue.put(tx)
                    logger.debug("Returning (%s) to Q. Q has (%s) items" % (tx, self.__main_queue.qsize()))
            self.__main_queue.task_done()

    @staticmethod
    def process_tx(data1, data2):
        if data1[0][0] != data2[0][0]:
            print("Error")
        else: tx=data1[0][0]
        logger.info("Processing %s" % (tx))
        return("Tx (%s) done. Len f1: (%s) - Len f2: (%s)" % (tx, len(data1), len(data2)))
        logger.info("Done processing %s" % (tx))
    
if __name__ == '__main__':
    d = events_combinator()
