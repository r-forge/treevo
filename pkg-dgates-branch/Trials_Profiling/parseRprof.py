#!/usr/bin/env python

import sys

#openfile= '/home/brightuser/Desktop/doRun.out'
#input1= open(openfile)
#S= input1.readlines()

class FuncCall():
    """Holds information on samples in a particular function, and the call hierarchy below it"""
    def __init__(self, name):
        self.name = name
        self.count = 1
        self.children = {}
        self.terminal_count = 0

    def incr(self):
        self.count += 1
    
    def add_children(self, child_list):
        """add children to the dictionary, and recursively add to those, etc."""
        if len(child_list) == 0:
            return
        try:
            if child_list[0] not in self.children:
                self.children[child_list[0]] = FuncCall(child_list[0])
            else:
                self.children[child_list[0]].incr()
            if len(child_list) > 1:
                self.children[child_list[0]].add_children(child_list[1:])
            else:
                self.children[child_list[0]].terminal_count += 1
        except:
            exit('problem adding child %s' % child_list[0])
        
    def output(self, lvl, tot):
        """write columnar output, with different levels of indentation for sucessively deeper func calls"""
        sys.stdout.write('%.5f%10d%10.5f%10d  ' % (self.count / float(tot), self.count, self.terminal_count / float(tot), self.terminal_count))
        for i in range(0, lvl):
            sys.stdout.write('| ')
        sys.stdout.write('%s\n' % (self.name))
        for child in sorted(list(self.children.values()), key=lambda c:c.count, reverse=True):
            child.output(lvl + 1, tot)
   

if len(sys.argv) != 2:
    exit('pass an output file from Rprof')

samples = [ line.strip().split() for line in open(sys.argv[1]) ]
tot_samples = len(samples)

toplevel = FuncCall('toplevel')
for samp in samples[1:]:
    samp.reverse()
    toplevel.add_children(samp)

sys.stdout.write('%7s%10s%10s%10s  %s\n' % ('tot%', 'tot#', 'self%', 'self#', 'function'))
#output= '/home/brightuser/Desktop/R_hierarchy'
#output=open(output,'w')
#output.write('%7s%10s%10s%10s  %s\n' % ('tot%', 'tot#', 'self%', 'self#', 'function'))
toplevel.output(0, tot_samples)

