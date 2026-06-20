
from collections import defaultdict


class ScalarCounter(object):
    def __init__(self, start=0):
        self.count = start


    def increment(self, step=1):
        self.count += step
        return self.count


    def decrement(self, step=1):
        return self.increment(step=-step)
            


class DictCounter(defaultdict):
    def __init__(self, start=dict(), step=1):
        defaultdict.__init__(self, int)
        for key in start:
            self.increment(key, start[key])

        
    def increment(self, key, step=1):
        self[key] += step
        return self[key]


    def decrement(self, key, step=1):
        return self.increment(key, -step)

    

class CountsTableRecord(list):
    def __init__(self, iterable=(), id=None):
        super().__init__(iterable)
        self.id = id

    @property
    def total(self):
        return sum(self)

    

class CountsTable(list):
    def __init__(self, iterable=(), samples=()):
        super().__init__(iterable)
        self.samples = samples

    def new_record(self, append=True):
        record = CountsTableRecord([0] * len(self.samples))
        if append:
            self.append(record)
        return record
