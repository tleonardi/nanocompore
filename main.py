import asyncio
import logging
from time import sleep
from aiofile import AIOFile, LineReader

logging.basicConfig(level=logging.DEBUG)

# setting up logger
logger = logging.getLogger("MyLogger")
logger.setLevel(logging.DEBUG)
console = logging.StreamHandler()
logger.addHandler(console)
MAX_CONNECTION = 1000


        
async def read_file(file, main_queue, data):
    logger.debug("Entered read_file with (%s)" % file)
    #for i in range(1000000000000):
    i=0
    for i in range(5):
        i+=1
        asyncio.sleep(10)
        res = {'Name':"TX(%s)"%i, 'Data':[["Name", "data", "data"]]}
        data[res['Name']] = res['Data']
        await main_queue.put(res["Name"])
        logger.debug("read_file: Added (%s) to Q. Q has (%s) items." % (res["Name"], main_queue.qsize()))
       
def process_tx(tx, data1, data2):
    data1.pop(tx)
    data2.pop(tx)
    print("processing (%s)"%tx)
    sleep(1)
    print("(%s) Done"%tx)
    results="Done"
    return results

async def consumer(main_queue, results, data1, data2):
    while True:
        logger.debug("Consumer - data1 has (%s) elements - data2 has (%s) elements." % (len(data1.keys()), len(data2.keys())))
        # Fetch a tx from the queue
        tx = await main_queue.get()
        logger.debug("consumer: picked (%s) from Q. Q has (%s) items." % (tx, main_queue.qsize()))
        # Check data is in data_1 and data_2, if not put back in main_queue

        if tx in data1.keys() and tx in data2.keys():
            results.append(await process_tx(tx, data1, data2))
            main_queue.task_done()
        else:
            await main_queue.put(tx)
            #logger.debug("consumer: Returning (%s) to Q. Q has (%s) items" % (tx, main_queue.qsize()))
            logger.debug(data1.keys())
            logger.debug(data2.keys())


async def consume_files(file1, file2, n=5):
    main_queue, results = asyncio.Queue(maxsize=20), []
    
    # we init the consumers, as the queues are empty at first,
    # they will be blocked on the main_queue.get()
    data_f1=dict()
    data_f2=dict()
    consumers = [asyncio.ensure_future(consumer(main_queue, results, data_f1, data_f2)) for _ in range(n)]

    producer1 = await read_file(file1, main_queue, data_f1)
    producer2 = await read_file(file2, main_queue, data_f2)
    #producers = [asyncio.ensure_future(read_file(file1, main_queue, data_f1)), read_file(file2, main_queue, data_f1)]
    logger.debug("Here")

    # wait for all item's inside the
    # main_queue to get task_done
    await main_queue.join()
    # cancel all coroutines
    #for consumer_future in consumers:
    #    consumer_future.cancel()
    return results


async def run():
    # we init more connectors to get better performance
    file1="/mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_kd/nanopolish/events.txt"
    file2="/mnt/home1/kouzarides/tl344/projects/nanopore_7SK/analysis_wt/nanopolish/events.txt"
    results = await consume_files(file1, file2, 100)
    print(results)


loop = asyncio.get_event_loop()
loop.run_until_complete(run())

