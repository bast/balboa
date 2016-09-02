"""
Unit tests.
"""
import sys
import os
import pytest


@pytest.fixture(scope='function')
def context(request):
    """
    Add context to test functions.
    """
    from aoeval import lib
    ctx = None
#   ctx = lib.aoeval_new()

#   def cleanup():
#       """
#       Clean up the context.
#       """
#       lib.aoeval_free(ctx)

#   request.addfinalizer(cleanup)
    return ctx


def test_aoeval(context):
    from aoeval import lib

    pass
