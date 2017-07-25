# -*- encoding: utf-8
"""
Tests for dummyforpy
"""
from numpy.testing import assert_array_equal
import numpy as np

from dummyforpy.core import calc_module


class TestSomething(object):
    """ A simple test """
    def test_mult_loop(self):
        f = 10
        a = np.arange(60).reshape(5, 4, 3)
        ref = a * f
        b = calc_module.calc(a, f)
        assert_array_equal(ref, b)
