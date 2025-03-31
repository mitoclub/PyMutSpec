import cProfile
import pstats
import sys


def profiler(_func=None, *, nlines=20):
    def decorator_profiler(func):
        def wrapper_profiler(*args, **kwargs):
            profile = cProfile.Profile()
            profile.enable()
            value = func(*args, **kwargs)
            profile.disable()
            ps = pstats.Stats(profile, stream=sys.stderr)
            print(f"{'-' * 30}\nFunction: {func.__name__}\n{'-' * 30}", file=sys.stderr)
            ps.sort_stats('cumtime', 'tottime')
            ps.print_stats(nlines)
            return value
        return wrapper_profiler

    if _func is None:
        return decorator_profiler
    else:
        return decorator_profiler(_func)
