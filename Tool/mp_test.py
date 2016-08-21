import time
from Tool import tot_multiprocess as tot_mp
from concurrent.futures import wait, ALL_COMPLETED


def mp_setup():
    mp = tot_mp.Multiprocessing(2)
    pool = mp.setup()
    return pool


def insignificant_task():
    time.sleep(10)
    return "I'm awake!"


def insignificant_callback(future):
    print(future.result())


def run():

    pool = mp_setup()

    futures = []  # Holding list for futures (jobs) submitted to process pool

    for x in range(3):
        job = pool.submit(insignificant_task)
        job.add_done_callback(insignificant_callback)
        futures.append(job)

    print("Start waiting: %s" % time.time())
    wait(futures, return_when=ALL_COMPLETED)
    print("Done waiting:  %s" % time.time())

    for x in range(3):
        job = pool.submit(insignificant_task)
        job.add_done_callback(insignificant_callback)
        futures.append(job)

    print("Second batch of jobs submitted!")

    print("About to tell Process Pool to close")
    pool.shutdown()  # Blocking
    print("Process Pool is now closed")


if __name__ == "__main__":
    run()
