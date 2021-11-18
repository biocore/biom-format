#!/usr/bin/env python

r"""Error profile

The error profile object to define different error types and handling within
BIOM. The following types are registered with the associated default states

empty : 'ignore'
    Treatment of an empty table (e.g., if Table.is_empty() is True). If a
    callable is provided, it implies 'call' and will set the callback function.

obssize : 'raise'
    Treatment of a table in which the number of observation ids does not match
    the size of the data.

sampsize : 'raise'
    Treatment of a table in which the number of sample ids does not match the
    size of the data.

obsdup : 'raise'
    Treatment of duplicate observation IDs.

sampdup : 'raise'
    Treatment of duplicate sample IDs.

obsmdsize : 'raise'
    Treatment of a table in which the number of observation metadata elements
    differs from the size of the data.

sampmdsize : 'raise'
    Treatment of a table in which the number of sample metadata elements
    differs from the size of the data.

Examples
--------

Use `seterr` to change how empty tables are handled:

>>> from numpy import inf
>>> from biom.err import seterr
>>> from biom import example_table
>>> make_empty_f = lambda v, i, md: v.sum() > inf
>>> _ = example_table.filter(make_empty_f, inplace=False)
>>> old_state = seterr(empty='raise')
>>> example_table.filter(make_empty_f, inplace=False)
Traceback (most recent call last):
...
TableException: Empty table!
>>> _ = seterr(**old_state)

"""

# -----------------------------------------------------------------------------
# Copyright (c) 2011-2020, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from warnings import warn
from sys import stdout
from contextlib import contextmanager

from biom.exception import TableException


# error messages
EMPTY = "Empty table!"
OBSSIZE = "Number of observation IDs differs from matrix size!"
SAMPSIZE = "Number of sample IDs differs from matrix size!"
OBSDUP = "Duplicate observation IDs"
SAMPDUP = "Duplicate sample IDs!"
OBSMDSIZE = "Size of observation metadata differs from matrix size!"
SAMPMDSIZE = "Size of sample metadata differs from matrix size!"


# _zz_ so the sort order places this test last
def _zz_test_empty(t):
    """Check if t is empty"""
    return t.is_empty()


def _test_obssize(t):
    """Check if the number of observations match data size"""
    return t.shape[0] != len(t.ids(axis='observation'))


def _test_sampsize(t):
    """Check if the number of samples match data size"""
    return t.shape[1] != len(t.ids(axis='sample'))


def _test_obsdup(t):
    """Check if there are duplicate observations"""
    return t.shape[0] != len(set(t.ids(axis='observation')))


def _test_sampdup(t):
    """Check if there are duplicate samples"""
    return t.shape[1] != len(set(t.ids(axis='sample')))


def _test_obsmdsize(t):
    """Check if the size of the observation metadata matches data size"""
    md = t.metadata(axis='observation')
    return t.shape[0] != len(md) if md is not None else False


def _test_sampmdsize(t):
    """Check if the size of the sample metadata matches data size"""
    md = t.metadata(axis='sample')
    return t.shape[1] != len(md) if md is not None else False


def _create_error_states(msg, callback, exception):
    """Create error states"""
    return {'ignore': lambda x: None,
            'warn': lambda x: warn(msg),
            'raise': lambda x: exception(msg),
            'call': callback if callback is not None else lambda x: None,
            'print': lambda x: stdout.write(msg + '\n')}


class ErrorProfile:
    """An error profile

    The error profile defines the types of errors that can be optionally
    handled, how those errors are handled, and performs the handling of the
    errors.
    """
    _valid_states = frozenset(['raise', 'ignore', 'call', 'print', 'warn'])

    def __init__(self):
        self._profile = {}
        self._state = {}
        self._test = {}

    def register(self, errtype, msg, state, test, callback=None,
                 exception=Exception):
        """Register an error type

        Paramters
        ---------
        errtype : str
            A name for the error, must be unique.
        msg : str
            An error message.
        state : {'raise', 'ignore', 'call', 'print', 'warn'}
            The default state.
        test : function
            A function to test for error.
        callback : function, optional
            A callback function for use with state 'call'
        exception : Exception, optional
            An exception to throw in state 'raises'.

        Raises
        ------
        KeyError
            If the errtype already exists
        KeyError
            If the state is invalid

        """
        if errtype in self:
            raise KeyError("Already registered: %s" % errtype)

        if state not in self._valid_states:
            raise KeyError("Unknown state: %s" % state)

        self._profile[errtype] = _create_error_states(msg, callback, exception)
        self._state[errtype] = state
        self._test[errtype] = test

    def unregister(self, errtype):
        """Unregister an error type

        Parameters
        ----------
        errtype : str
            The error type to unregister

        Raises
        ------
        KeyError
            If the error type is unknown

        Returns
        -------
        dict
            The profile associated with the error.
        function
            The test function associated with the error.
        str
            The state at time of unregistering associated with the error

        """
        if errtype not in self:
            raise KeyError("Unknown error type: %s" % errtype)

        prof = self._profile.pop(errtype)
        func = self._test.pop(errtype)
        state = self._state.pop(errtype)

        return (prof, func, state)

    @property
    def state(self):
        """Return current state"""
        return self._state

    @state.setter
    def state(self, new_state):
        """Update current state"""
        if 'all' in new_state:
            to_update = [(err, new_state['all']) for err in self._state]
        else:
            to_update = new_state.items()

        for errtype, new_state in to_update:
            if new_state not in self._valid_states:
                raise KeyError("Unknown state type: %s" % new_state)
            if errtype not in self._state:
                raise KeyError("Unknown error type: %s" % errtype)

            self._state[errtype] = new_state

    def __contains__(self, errtype):
        """Check if an error type exists"""
        return errtype in self._state

    def test(self, item, *args):
        """Test for an error

        Parameters
        ----------
        item : object
            An item to test for error
        *args : list, optional
            Error types to check, if not provided, all known error types are
            checked.

        Examples
        --------
        >>> from biom import example_table
        >>> from biom.err import __errprof
        >>> __errprof.test(example_table, 'empty')

        """
        if not args:
            args = self._test.keys()

        for errtype in sorted(args):
            test = self._test.get(errtype, lambda: None)

            if test(item):
                return self._handle_error(errtype, item)

    def _handle_error(self, errtype, item):
        """Handle an error"""
        state = self._state[errtype]
        profile = self._profile[errtype]
        return profile[state](item)

    def setcall(self, errtype, func):
        """Set a function callback

        Parameters
        ----------
        errtype : str
            The error type to set a callback for.
        func : function
            The callback function. This function must accept a single
            parameter. The result of this function will be returned through
            the handler.

        Raises
        ------
        KeyError
            If an unknown error type is specified.

        Examples
        --------
        >>> from biom.err import __errprof
        >>> oldcall = __errprof.setcall('empty', lambda item: 123)
        >>> _ = __errprof.setcall('empty', oldcall)

        """
        if errtype not in self:
            raise KeyError("Unknown error type: %s" % errtype)

        old_call = self._profile[errtype]['call']
        self._profile[errtype]['call'] = func
        return old_call

    def getcall(self, errtype):
        """Get a function callback

        Parameters
        ----------
        errtype : str
            The error type to get the callback for.

        Raises
        ------
        KeyError
            If the error type is unknown.

        Returns
        -------
        function
            The associated callback function.

        """
        if errtype not in self:
            raise KeyError("Unknown error type: %s" % errtype)

        return self._profile[errtype]['call']


__errprof = ErrorProfile()
__errprof.register('empty', EMPTY, 'ignore', _zz_test_empty,
                   exception=TableException)
__errprof.register('obssize', OBSSIZE, 'raise', _test_obssize,
                   exception=TableException)
__errprof.register('sampsize', SAMPSIZE, 'raise', _test_sampsize,
                   exception=TableException)
__errprof.register('obsdup', OBSDUP, 'raise', _test_obsdup,
                   exception=TableException)
__errprof.register('sampdup', SAMPDUP, 'raise', _test_sampdup,
                   exception=TableException)
__errprof.register('obsmdsize', OBSMDSIZE, 'raise', _test_obsmdsize,
                   exception=TableException)
__errprof.register('sampmdsize', SAMPMDSIZE, 'raise', _test_sampmdsize,
                   exception=TableException)


def geterr():
    """Returns the current error profile state"""
    return __errprof.state.copy()


def seterr(**kwargs):
    """How table errors are handled, API based on numpy's seterr

    Notes
    ----------
    all : {'ignore', 'warn', 'raise', 'call', 'print'}, optional
        Set treatment for all error types

        - ignore: Take no action when the exception occurs.
        - warn: Print a `RuntimeWarning` (via the Python `warnings` module).
        - raise: Raise a `TableException`.
        - call: Call a function specified using the `seterrcall` function.
        - print: Print a warning directly to ``stdout``.

    Returns
    -------
    old_settings: dict
        Dictionary containing the old settings

    See also
    --------
    geterr
    seterrcall
    geterrcall
    errstate

    Examples
    --------

    >>> from numpy import inf
    >>> from biom.err import seterr
    >>> from biom import example_table
    >>> make_empty_f = lambda v, i, md: v.sum() > inf
    >>> _ = example_table.filter(make_empty_f, inplace=False)
    >>> old_state = seterr(empty='raise')
    >>> example_table.filter(make_empty_f, inplace=False)
    Traceback (most recent call last):
    ...
    TableException: Empty table!
    >>> _ = seterr(**old_state)

    """
    old_state = __errprof.state.copy()
    if 'all' in kwargs:
        __errprof.state = {'all': kwargs['all']}
    else:
        __errprof.state = kwargs
    return old_state


def seterrcall(errtype, func):
    """Set an error callback function

    This is similar to numpy's seterrcall, except that we define it by error
    type. The justification is that numpy's errors are focused specifically on
    floating point errors, while the types of errors we support in BIOM are
    more broad.

    Parameters
    ----------
    errtype : str
        Must be a valid key in biom.__errprof
    func : callable
        A function that is called if the given error is encountered and if the
        current error profile indicates 'call' for the error type.

    Returns
    -------
    callable or None
        The old error handler

    Raises
    ------
    KeyError
        If the errtype specified is not in the error profile.

    """
    if errtype not in __errprof:
        raise KeyError("Unknown error type: %s" % errtype)
    else:
        return __errprof.setcall(errtype, func)


def geterrcall(errtype):
    """Get the current callback function associated with an error type

    Parameters
    ----------
    errtype : str
        Must be a valid key in the error profile

    Returns
    -------
    callable
        The current callback function

    Raises
    ------
    KeyError
        If the errtype specified is not in the error profile.

    """
    if errtype not in __errprof:
        raise KeyError("Unknown error type: %s" % errtype)
    else:
        return __errprof.getcall(errtype)


def errcheck(table, *errtypes):
    """Check if there is an error, and respond appropriately

    Parameters
    ----------
    table : Table
        The table to check
    errtypes : vargs of str, if not specified, defaults to all errors
        Errors to test

    Notes
    -----
    If an error type has 'raises' as its state, the exception is raised from
    here to avoid cluttering the traceback.

    Returns
    -------
    dependent on error state and setting

    """
    ret = __errprof.test(table, *errtypes)
    if isinstance(ret, Exception):
        raise ret
    else:
        return ret


@contextmanager
def errstate(**kwargs):
    """Context manager for error handling

    Using an instance of `errstate` as a context manager allows statements in
    that context to execute with a known error handling behavior. Upon entering
    the context the error handling is set with `seterr`, and upon exiting it is
    reset to what it was before. Please note, this text was taken near verbatim
    from numpy's errstate method.

    Parameters
    ----------
    kwargs : {empty}
        Keyword arguments. The valid error types that are defined. Each keyword
        should have a string or callable for the particular error. Values are:
        {'ignore', 'warn', 'raise', 'call', 'print'}

    See Also
    --------
    seterr, geterr, seterrcall, geterrcall

    """
    old_state = seterr(**kwargs)
    yield
    seterr(**old_state)
