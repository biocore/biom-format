#!/usr/bin/env python

import numpy
from numpy import testing, array, zeros, asarray
from unittest import main, TestSuite, findTestCases, \
    TestCase as TestCaseOriginal

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, BIOM-Format Project"
__credits__ = ["Rob Knight", "Peter Maxwell", "Sandra Smit",
               "Zongzhi Liu", "Micah Hamady", "Daniel McDonald",
               "Jai Ram Rideout"]
__license__ = "GPL"
__url__ = "http://biom-format.org"
__version__ = "1.1.2"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"

class TestCase(TestCaseOriginal):
    """Add in some support for numpy vectors

    Methods pulled from PyCogent (http://pycogent.sourceforge.net)
    """
    def assertFloatEqualRel(self, obs, exp, eps=1e-6):
        """Tests whether two floating point numbers/arrays are approx. equal.

        Checks whether the distance is within epsilon relative to the value
        of the sum of observed and expected. Use this method when you expect
        the difference to be small relative to the magnitudes of the observed
        and expected values.

        Note: for arbitrary objects, need to compare the specific attribute
        that's numeric, not the whole object, using this method.
        """
        #do array check first
        #note that we can't use array ops to combine, because we need to check
        #at each element whether the expected is zero to do the test to avoid
        #floating point error.
        #WARNING: numpy iterates over objects that are not regular Python
        #floats/ints, so need to explicitly catch scalar values and prevent
        #cast to array if we want the exact object to print out correctly.
        is_array = False
        if hasattr(obs, 'keys') and hasattr(exp, 'keys'):   #both dicts?
            result = self._get_values_from_matching_dicts(obs, exp)
            if result:
                obs, exp = result
        else:
            try:
                iter(obs)
                iter(exp)
            except TypeError:
                obs = [obs]
                exp = [exp]
            else:
                try:
                    arr_obs = array(obs)
                    arr_exp = array(exp)
                    arr_diff = arr_obs - arr_exp
                    if arr_obs.shape != arr_exp.shape:
                        self.fail("Wrong shape: Got %s, but expected %s" % \
                            (`obs`, `exp`))
                    obs = arr_obs.ravel()
                    exp = arr_exp.ravel()
                    is_array=True
                except (TypeError, ValueError):
                    pass

        # shape mismatch can still get by...
        # explict cast is to work around bug in certain versions of numpy
        # installed version on osx 10.5
        if asarray(obs, object).shape != asarray(exp, object).shape:
            self.fail("Wrong shape: Got %s, but expected %s" % (obs, exp))

        for observed, expected in zip(obs, exp):
            #try the cheap comparison first
            if observed == expected:
                continue
            try:
                sum = float(observed + expected)
                diff = float(observed - expected)
                if (sum == 0):
                    if is_array:
                        self.failIf(abs(diff) > abs(eps), \
                            "Got %s, but expected %s (diff was %s)" % \
                            (`arr_obs`, `arr_exp`, `arr_diff`))
                    else:
                        self.failIf(abs(diff) > abs(eps), \
                            "Got %s, but expected %s (diff was %s)" % \
                            (`observed`, `expected`, `diff`))

                else:
                    if is_array:
                        self.failIf(abs(diff/sum) > abs(eps), \
                            "Got %s, but expected %s (diff was %s)" % \
                            (`arr_obs`, `arr_exp`, `arr_diff`))
                    else:
                        self.failIf(abs(diff/sum) > abs(eps), \
                            "Got %s, but expected %s (diff was %s)" % \
                            (`observed`, `expected`, `diff`))
            except (TypeError, ValueError, AttributeError, NotImplementedError):
                self.fail("Got %s, but expected %s" % \
                    (`observed`, `expected`))

    def assertFloatEqualAbs(self, obs, exp, eps=1e-6):
        """
        Tests whether two floating point numbers are approximately equal.

        Checks whether the absolute value of (a - b) is within epsilon. Use
        this method when you expect that one of the values should be very
        small, and the other should be zero.
        """
        #do array check first
        #note that we can't use array ops to combine, because we need to check
        #at each element whether the expected is zero to do the test to avoid
        #floating point error.
        if hasattr(obs, 'keys') and hasattr(exp, 'keys'):   #both dicts?
            result = self._get_values_from_matching_dicts(obs, exp)
            if result:
                obs, exp = result
        else:
            try:
                iter(obs)
                iter(exp)
            except TypeError:
                obs = [obs]
                exp = [exp]
            else:
                try:
                    arr_obs = array(obs)
                    arr_exp = array(exp)
                    if arr_obs.shape != arr_exp.shape:
                        self.fail("Wrong shape: Got %s, but expected %s" % \
                            (`obs`, `exp`))
                    diff = arr_obs - arr_exp
                    self.failIf(abs(diff).max() > eps, \
                        "Got %s, but expected %s (diff was %s)" % \
                        (`obs`, `exp`, `diff`))
                    return
                except (TypeError, ValueError):
                    pass
        #only get here if array comparison failed
        for observed, expected in zip(obs, exp):
            #cheap comparison first
            if observed == expected:
                continue
            try:
                diff = observed - expected
                self.failIf(abs(diff) > abs(eps),
                        "Got %s, but expected %s (diff was %s)" % \
                        (`observed`, `expected`, `diff`))
            except (TypeError, ValueError, AttributeError, NotImplementedError):
                self.fail("Got %s, but expected %s" % \
                    (`observed`, `expected`))

    def assertFloatEqual(self, obs, exp, eps=1e-6, rel_eps=None, \
                         abs_eps=None):
        """Tests whether two floating point numbers are approximately equal.

        If one of the arguments is zero, tests the absolute magnitude of the
        difference; otherwise, tests the relative magnitude.

        Use this method as a reasonable default.
        """
        obs = numpy.asarray(obs, dtype='O')
        exp = numpy.asarray(exp, dtype='O')
        obs = numpy.ravel(obs)
        exp = numpy.ravel(exp)

        if obs.shape != exp.shape:
            self.fail("Shape mismatch. Got, %s but expected %s" % (obs, exp))

        for observed, expected in zip(obs, exp):
            if self._is_equal(observed, expected):
                continue
            try:
                rel_eps = rel_eps or eps
                abs_eps = abs_eps or eps
                if (observed == 0) or (expected == 0):
                    self.assertFloatEqualAbs(observed, expected, abs_eps)
                else:
                    self.assertFloatEqualRel(observed, expected, rel_eps)
            except (TypeError, ValueError, AttributeError, NotImplementedError):
                self.fail("Got %s, but expected %s" % \
                        (`observed`, `expected`))


    def _is_equal(self, observed, expected):
        """Returns True if observed and expected are equal, False otherwise."""
        #errors to catch: TypeError when obs is None
        tolist_errors = (AttributeError, ValueError, TypeError)

        try:
            obs = observed.tolist()
        except tolist_errors:
            obs = observed
        try:
            exp = expected.tolist()
        except tolist_errors:
            exp = expected
        return obs == exp

    def failUnlessEqual(self, observed, expected, msg=None):
        """Fail if the two objects are unequal as determined by !=

        Overridden to make error message enforce order of observed, expected.
        Use numpy.testing.assert_equal if ValueError, TypeError raised.
        """
        try:
            if not self._is_equal(observed, expected):
                raise self.failureException, \
                (msg or 'Got %s, but expected %s' % (`observed`, `expected`))
        except (ValueError, TypeError), e:
            #The truth value of an array with more than one element is
            #ambiguous. Use a.any() or a.all()
            #descriptor 'tolist' of 'numpy.generic' object needs an argument
            testing.assert_equal(observed, expected)
    #following needed to get our version instead of unittest's
    assertEqual = assertEquals = failUnlessEqual

    def assertNotSameObj(self, observed, expected, msg=None):
        """Fail if 'observed is expected'"""
        try:
            if observed is not expected:
                return
        except:
            pass
        raise self.failureException, \
        (msg or 'Observed %s is the same as expected %s' % \
        (`observed`, `expected`))
