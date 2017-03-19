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
        start = time.clock()
        try:
            return f(*args, **kwargs)
        except:
            raise
        finally:
            end = time.clock()
            key = f.__name__
            counted.called[key] += 1
            counted.timing[key] += end - start
    return wrapper

counted.called = Counter()
counted.timing = Counter()
