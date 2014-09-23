#!/usr/bin/env python

"""
havaGen is for generating sequence files of the imaginary type 'hava'.
These files consist of attributes in the following newline seperated order
hava
quarzk
muchalo
fakours
selimizicka
marshoon

Since this 'sequence' has absolutely no utility outside of screed, it's only
purpose is to make sure screed can work with arbitrary fields when running
the nosetests.

This is a work of fiction. Names are the product of the author's imagination
and any resemblance to real life is entirely coincidental.
"""

import sys, os
import random

class collectionOFiles(object):
    def __init__(self, baseName, divisions, totalSize):
        self.baseName = baseName
        self.divisions = divisions
        self.totalSize = totalSize

        self.fileHandles = {}
        for i in range(0, divisions):
            filename = self.baseName + "_%d" % i
            fh = open(filename, "wb")
            divisor = i * 2
            if divisor == 0:
                divisor = 1
            self.fileHandles[filename]= (fh, self.totalSize/divisor, 0)

    def writeRecord(self, hava, quarzk, muchalo, fakours, selimizicka, marshoon):
        toRemove = []
        for filename in self.fileHandles:
            file, limit, count = self.fileHandles[filename]
            file.write("%s\n%s\n%s\n%s\n%s\n%s\n" % (hava, quarzk, muchalo, fakours, selimizicka, marshoon))
            count += 1
            if count >= limit:
                file.close()
                toRemove.append(filename)
            else:
                self.fileHandles[filename] = (file, limit, count)

        for fh in toRemove:
            self.fileHandles.pop(fh)

    def finished(self):
        return len(self.fileHandles) == 0

def genString(length, allowedChars):
    res = []
    for i in range(0, length):
        char = allowedChars[random.randint(0, len(allowedChars)-1)]
        res.append(char)
    return "".join(res)

def createHavaFiles(filename, size, divisions):
    cof = collectionOFiles(filename, divisions, size)
    counter = 0
    lenString = 80
    allowedQuarzk = ['A', 'T', 'C', 'G']
    allowedMuchalo = "A B C D E F G H I J K L M N O P".split(' ')
    allowedFakours = "1 2 3 4 5 6 7 8 9".split(' ')
    allowedSelimizicka = ["b"]
    allowedMarshoon = "A 1 B 2 C 3 D 4 E 5 G 6 F 7".split(' ')
    while(not cof.finished()):
        hava = "test_00%d" % counter
        quarzk = genString(lenString, allowedQuarzk)
        muchalo = genString(lenString, allowedMuchalo)
        fakours = genString(lenString, allowedFakours)
        selimizicka = genString(lenString, allowedSelimizicka)
        marshoon = genString(lenString, allowedMarshoon)
        cof.writeRecord(hava, quarzk, muchalo, fakours, selimizicka, marshoon)
        counter += 1
    return

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print "Usage: <filename> <size> <divisions>"
        exit(1)

    filename = sys.argv[1]
    size = int(sys.argv[2])
    divisions = int(sys.argv[3])

    createHavaFiles(filename, size, divisions)
