"""Minimal stand-in for pytest, so the repository's own fixtures and report
builders can be driven outside a test run. The fixture functions in
tests/conftest.py and tests/test_merge_report.py are plain functions once the
decorator is removed; this module removes it.

Only what those two modules import is provided.
"""


def fixture(*args, **kwargs):
    if args and callable(args[0]) and not kwargs:
        return args[0]

    def deco(fn):
        return fn
    return deco


def skip(reason=""):
    raise RuntimeError("skip: %s" % reason)


def fail(reason=""):
    raise AssertionError(reason)


class _Mark:
    def parametrize(self, *a, **k):
        def deco(fn):
            return fn
        return deco

    def __getattr__(self, _name):
        def deco(fn=None, *a, **k):
            return fn if fn is not None else (lambda f: f)
        return deco


mark = _Mark()
