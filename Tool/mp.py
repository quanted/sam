from concurrent.futures import ProcessPoolExecutor as Pool
from concurrent.futures import wait
from concurrent.futures import ALL_COMPLETED
# import concurrent.futures


class Multiprocessing:
    def __init__(self, nproc):
        self.nproc = nproc

    def setup(self):
        if not self.nproc:
            import multiprocessing
            self.nproc = multiprocessing.cpu_count()

        print("CPU cores = %s" % self.nproc)

        return Pool(max_workers=int(self.nproc))

    @staticmethod
    def wait_to_finish(futures):
        wait(futures, return_when=ALL_COMPLETED)

    def chunk(self, no_upstream_reaches):
        chunk_size = int(no_upstream_reaches / self.nproc)

        # Generate a list of evenly distributed ranges to chunk the list of reaches by
        return self.reach_ranges(chunk_size, no_upstream_reaches)

    def reach_ranges(self, chunk_size, no_upstream_reaches):
        """
        Creates a list of ranges for the 'no_upstream_reaches' with a length equal to the 'nproc'.  Each range represents
        a sequential section of 'no_upstream_reaches' with a length of 'chunk_size'.

        :param chunk_size:
        :param nproc:
        :param no_upstream_reaches:
        :return: list of tuples, len = nproc, tuple = (range_0, range_n)
        """
        reach_ranges = []

        for x in range(self.nproc):
            if x == 0:
                reach_ranges.append((0, chunk_size))
            elif x < self.nproc - 1:
                reach_ranges.append((chunk_size * x, (chunk_size * x) + chunk_size))
            else:
                reach_ranges.append(
                    (chunk_size * x,
                     (chunk_size * x) + (no_upstream_reaches - (chunk_size * (self.nproc - 1))))
                )

        return reach_ranges

    # @staticmethod
    # def reach_calc(reaches):
    #     """
    #     Function to send to each process when using multiprocessing. The number total number of reaches processed is
    #     first sent through reach_ranges_mp() to generate ranges to divide up the work that is sent into this function.
    #     This function is a batch of reaches to be processed in a single process.
    #     :param reaches: list of Reach class instances
    #     """
    #     for reach in reaches:
    #         if not reach.run:
    #             reach.process()
