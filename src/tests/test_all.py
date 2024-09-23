import pytest
import mavdac


def test_sum_as_string():
    assert mavdac.sum_as_string(1, 1) == "2"
