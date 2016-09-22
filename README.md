[![Build Status](https://travis-ci.org/bast/balboa.svg?branch=master)](https://travis-ci.org/bast/balboa/builds)

# balboa

Experimental code. Under heavy development, nothing is stable.

### Testing

```
PROJECT_BUILD_DIR=$PWD/build PROJECT_INCLUDE_DIR=$PWD/api PYTHONPATH=$PWD/api py.test -vv test/test.py
```

### TODO

- Get rid of AO_BLOCK_LENGTH; make block length flexible
- Describe installation, configuration, and API

### Ordering of AOs

Let us assume that we have N basis functions and P points.


```
[ geo_000 (undifferentiated)         ][ geo_100 ][ geo_010 ][ geo_001 ] ...
[ ao_1    ][ ao_2    ] ... [ ao_3    ]
[ 1 ... P ][ 1 ... P ] ... [ 1 ... P ]
```
