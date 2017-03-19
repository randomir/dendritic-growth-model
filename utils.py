from functools import wraps


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
