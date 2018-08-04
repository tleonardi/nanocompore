import logging
from time import sleep
from queue import Queue
import threading
import mmap

logging.basicConfig(level=logging.DEBUG)

# setting up logger
logger = logging.getLogger("Logger")
logger.setLevel(logging.DEBUG)
console = logging.StreamHandler()


def mmap_readline(file, sep="\t"):
    with open(file, "r+b") as f:
        m = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
        while True:
            line=m.readline()
            if line == b'': break
            yield line.decode().split(sep)[0:7]

def read_file2(file, main_queue, data):
    tx="contig"
    res=[]
    for line in mmap_readline(file):
        if line[0]==tx:
            res.append(line)
            tx=line[0]
            logger.debug("Adding line to (%s)" % (tx))
        else:
            data[tx]=res
            main_queue.put(tx)
            res=line
            tx=line[0]
            logger.debug("Finished processing (%s), (%s) lines" % (tx, len(data[tx])))
       
def read_file(file, main_queue, data):
    logger.debug("Entered read_file with (%s)" % file)
    #for i in range(1000000000000):
    i=0
    for i in range(5):
        i+=1
        sleep(0.3)
        res = {'Name':"TX(%s)"%i, 'Data':[["Name", "data", "data"]]}
        data[res['Name']] = res['Data']
        main_queue.put(res["Name"])
        logger.debug("read_file: Added (%s) to Q. Q has (%s) items." % (res["Name"], main_queue.qsize()))
       
def process_tx(tx, data1, data2):
    data1.pop(tx)
    data2.pop(tx)
    print("processing (%s)"%tx)
    sleep(10)
    print("(%s) Done"%tx)
    results="Done"
    return results

def consumer(main_queue, results, data1, data2, processed):
    while True:
        # Fetch a tx from the queue
        tx =  main_queue.get()
        if tx is None:
            break
        logger.debug("consumer: picked (%s) from Q. Q has (%s) items." % (tx, main_queue.qsize()))
        # Check data is in data_1 and data_2, if not put back in main_queue
        if tx not in processed:
            if tx in data1.keys() and tx in data2.keys():
                results.append( process_tx(tx, data1, data2))
                processed.add(tx)
            else:
                main_queue.put(tx)
                logger.debug("consumer: Returning (%s) to Q. Q has (%s) items" % (tx, main_queue.qsize()))
        main_queue.task_done()



main_queue = Queue(maxsize=20)
results=[] 
data_f1=dict()
data_f2=dict()

file1="/mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_wt/nanopolish/events.txt"
file2="/mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_wt/nanopolish/events.txt"

processed=set()
reader1 = threading.Thread(name="reader1", target=read_file2, args=(file1, main_queue, data_f1,))
#reader2 = threading.Thread(name="reader2", target=read_file2, args=(file2, main_queue, data_f2,))

consumers=[]
for i in range(3):
    consumers.append(threading.Thread(name="consumer (%s)"%i, target=consumer, args=(main_queue, results, data_f1, data_f2, processed)))

reader1.start()
#reader2.start()

#for c in consumers:
    #c.start()

reader1.join()
#reader2.join()
print(data_f1)
for _ in consumers:
    main_queue.put(None)
