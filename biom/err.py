from warnings import warn
from sys import stdout
from contextlib import contextmanager
import types

from biom.exception import TableException


def _raise_ex(ex):
    """Raise an exception, useful for raising from lambdas"""
    raise ex


class _ErrorProfile(object):
    _valid_states = frozenset(['raise', 'ignore', 'call', 'print', 'warn'])

    _profile = {'emptytable':
                {'ignore': lambda: None,
                 'warn': lambda: warn('Empty table!'),
                 'raise': lambda: _raise_ex(TableException),
                 'call': lambda: None,
                 'print': lambda: stdout.write('Empty table!')}}

    _state = {'emptytable': 'raise'}
    _test = (('emptytable', lambda t: t.isempty()),)

    @property
    def state(self):
        """Return current state"""
        return self._state

    @state.setter
    def state(self, **kwargs):
        """Update current state"""
        for errtype, new_state in kwargs.items():
            if new_state not in self._valid_states:
                raise KeyError("Unknown state type: %s" % new_state)
            self._state[errtype] = new_state

    def __getitem__(self, errtype):
        """Get the current profile for an errortype"""
        return self._profile[errtype]

    def __contains__(self, errtype):
        """Check if an error type exists"""
        return errtype in self._state

    def test(self, item):
        """Test an item for error"""
        for errtype, test_func in self._test:
            if self._state[errtype] != 'ignore' and test_func(item):
                self._handle_error(errtype)

    def _handle_error(self, errtype):
        """Handle an error"""
        state = self._state[errtype]
        self._profile[errtype][state]()

    def setcall(self, errtype, func):
        """Set a function callback"""
        old_call = self._profile[errtype]['call']
        self._profile[errtype]['call'] = func
        return old_call

    def getcall(self, errtype):
        """Get a function callback"""
        return self._profile[errtype]['call']


__error_profile = _ErrorProfile()


def geterr():
    """Returns the current error profile state"""
    return __error_profile.state


def seterr(all=None, emptytable=None):
    """How table errors are handled, API based on numpy's seterr

    Parameters
    ----------
    all : {'ignore', 'warn', 'raise', 'call', 'print', callable}, optional
        Set treatment for all error types

        - ignore: Take no action when the exception occurs.
        - warn: Print a `RuntimeWarning` (via the Python `warnings` module).
        - raise: Raise a `TableException`.
        - call: Call a function specified using the `seterrcall` function.
        - print: Print a warning directly to ``stdout``.
        - callable: implies 'call' and will set the callback function for all
            error types.

    emptytable : {'ignore', 'warn', 'raise', 'call', 'print', callable},
        optional
        Treatment of an empty table (e.g., Table.sum(axis='whole') == 0). If a
        callable is provided, it implies 'call' and will set the callback
        function.

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
    >>> from biom.util import seterr
    >>> from biom import example_table
    >>> example_table.filter(lambda v, i, md: v.sum() > inf, inplace=False)
    Traceback (most recent call last):
    ...
    biom.exception.TableException: All data was filtered out!
    >>> seterr(all='ignore')
    >>> example_table.filter(lambda v, i, md: v.sum() > inf, inplace=False)

    """
    old_state = __error_profile.state.copy()
    if all is not None:
        for k in __error_profile.state:
            __error_profile.state[k] = all
    elif emptytable is not None:
        if isinstance(emptytable, types.FunctionType):
            __error_profile.setcall('emptytable', emptytable)
        else:
            __error_profile.state['emptytable'] = emptytable
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
        Must be a valid key in biom.__error_profile
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
    if errtype not in __error_profile:
        raise KeyError("Unknown error type: %s" % errtype)
    else:
        return __error_profile.setcall(errtype, func)


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
    if errtype not in __error_profile:
        raise KeyError("Unknown error type: %s" % errtype)
    else:
        return __error_profile.getcall(errtype)


def check_error(table):
    """Check if there is an error, and respond appropriately

    Parameters
    ----------
    table : Table
        The table to check

    Returns
    -------
    dependent on error state and setting

    """
    __error_profile.test(table)


@contextmanager
def errstate(**kwargs):
    """Context manager for error handling

    Using an instance of `errstate` as a context manager allows statements in
    that context to execute with a known error handling behavior. Upon entering
    the context the error handling is set with `seterr` and `seterrcall`, and
    upon exiting it is reset to what it was before. Please note, this text was
    taken verbatim from numpy's errstate method.

    Parameters
    ----------
    kwargs : {emptytable}
        Keyword arguments. The valid error types that are defined. Each keyword
        should have a string or callable for the particular error. If the value
        is a callable, it implies a callback is to be specified for the error
        type. Values are: {'ignore', 'warn', 'raise', 'call', 'print',
                           callable}

    See Also
    --------
    seterr, geterr, seterrcall, geterrcall

    """
    old_state = seterr(**kwargs)
    yield
    seterr(old_state)
