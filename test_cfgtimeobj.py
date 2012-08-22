#!/usr/bin/env python
""" Test suite"""

import unittest
import configtimeobj
import math

class TestCorrelator(unittest.TestCase):

    def setUp(self):
        self.example = {'a': {0:1.0, 1:2.0}, 'b': {0:3.0, 1:4.0 } , 'c': {0:5.0, 1:0.0 }}
        self.cto = configtimeobj.Cfgtimeobj.fromDataDict(self.example)
        
    def test_setup(self):
        self.assertEqual(self.example['a'],self.cto.get(config='a'))
        self.assertEqual(self.example['b'],self.cto.get(config='b'))        
        self.assertEqual(self.example['b'][0],self.cto.get(config='b',time=0))                

    def test_bad_setup(self):
        self.assertRaises(Exception,configtimeobj.Cfgtimeobj.fromDataDict,({'a': {0:1.0 }, 'b': {0:5.0, 1:4.0 }}))
        self.assertRaises(Exception,configtimeobj.Cfgtimeobj.fromDataDict,({'a': {0:1.0, 2:2.0 }, 'b': {0:5.0, 1:4.0 }}))


    def test_averages(self):
        self.assertEqual(self.cto.average_over_times(),{'a': 1.5, 'b': 3.5, 'c':2.5 })
        self.assertEqual(self.cto.sum_over_configs(), {0: 9.0, 1: 6.0 })
        self.assertEqual(self.cto.average_over_configs(), {0: 3.0, 1: 2.0 })
        self.assertEqual(self.cto.average_all(), 2.5)        
        self.assertEqual(self.cto.jackknifed_averages()['a'], {0: 4.0, 1: 2.0 })
        self.assertEqual(self.cto.jackknifed_averages()['b'], {0: 3.0, 1: 1.0 })    
        self.assertEqual(self.cto.jackknifed_averages()['c'], {0: 2.0, 1: 3.0 })

        self.assertEqual(self.cto.jackknifed_full_average()['a'], 3.0)
        self.assertEqual(self.cto.jackknifed_full_average()['b'], 2.0)
        self.assertEqual(self.cto.jackknifed_full_average()['c'], 2.5)
    # def test_errorbars(self):
    #     self.assertEqual(self.cto.jackknifed_errors(), {0: math.sqrt(float(4)/float(3)), 1:  math.sqrt(float(4)/float(3)) })

    # def test_average_sub_vev(self):
    #     self.assertEqual(self.cto.average_vev(), (1.5+3.5+2.5)/3.0) # 2.5
    #     self.assertEqual(self.cto.average_sub_vev(), {0:0.5, 1:-0.5}  )
        
if __name__ == '__main__':
    unittest.main()
