from heapq import heappush, heappop

class RunningMedian:
    def __init__(self):
        self.lowers, self.highers = [], []

    def add_number(self, number):
        if not self.highers or number > self.highers[0]:
            heappush(self.highers, number)
        else:
            heappush(self.lowers, -number)  # for lowers we need a max heap
        self.rebalance()

    def rebalance(self):
        if len(self.lowers) - len(self.highers) > 1:
            heappush(self.highers, -heappop(self.lowers))
        elif len(self.highers) - len(self.lowers) > 1:
            heappush(self.lowers, -heappop(self.highers))

    def get_median(self):
        if len(self.lowers) == len(self.highers):
            return (-self.lowers[0] + self.highers[0])/2
        elif len(self.lowers) > len(self.highers):
            return -self.lowers[0]
        else:
            return self.highers[0]