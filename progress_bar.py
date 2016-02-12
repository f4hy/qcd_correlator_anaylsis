#!/usr/bin/env python2
import sys

class progress_bar(object):

    def __init__(self, total=100):

        self.total = total

    def update(self, progress):
        percent = 100*progress/self.total
        if self.total > 1:
            sys.stdout.write('\r[{0}] {1}/{2}'.format('#'*(percent/10)+' '*(10-percent/10), progress, self.total))
            sys.stdout.flush()

    def done(self):
        if self.total > 1:
            sys.stdout.write('\r[{0}] {1}/{2} DONE!\n'.format('#'*10, self.total, self.total))


def test():
    import time

    p = progress_bar()

    # p.update(10)
    # p.update(20)

    for i in range(0,100,10):
        time.sleep(1)
        p.update(i)
    p.done()

if __name__ == "__main__":
    test()
