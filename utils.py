import time
from functools import wraps
from collections import Counter


def memoized(f):
    """Function memoization. Both positional and keyword arguments are supported."""
    memo = {}
    @wraps(f)
    def wrapper(*args, **kwargs):
        key = repr((args, kwargs))
        if not key in memo:
            memo[key] = f(*args, **kwargs)
        return memo[key]
    wrapper.uncache = memo.clear
    return wrapper


def counted(f):
    """Count function calls. Simple profiling."""
    @wraps(f)
    def wrapper(*args, **kwargs):
        # for python>3.3, use ``time.perf_counter``
        start = time.clock()
        try:
            return f(*args, **kwargs)
        except Exception as e:
            counted.exceptions[repr(e)] += 1
            raise
        finally:
            end = time.clock()
            counted.called[f.__name__] += 1
            counted.timing[f.__name__] += end - start
    return wrapper

counted.called = Counter()
counted.timing = Counter()
counted.exceptions = Counter()
