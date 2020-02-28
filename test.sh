#!/usr/bin/env bash

export BALBOA_LIBRARY_DIR=$PWD/build/lib
export BALBOA_INCLUDE_DIR=$PWD/balboa
export PYTHONPATH=$PWD

pytest -vv -s test/test.py
