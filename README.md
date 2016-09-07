[![Build Status](https://travis-ci.org/bast/balboa.svg?branch=master)](https://travis-ci.org/bast/balboa/builds)

# balboa

Experimental code. Under heavy development, nothing is stable.

### Testing

```
PROJECT_BUILD_DIR=$PWD/build PROJECT_INCLUDE_DIR=$PWD/api PYTHONPATH=$PWD/api py.test -vv test/test.py
```

### TODO

- Max angular momentum and max geo diff order should be set by CMake
- Get rid of AO_BLOCK_LENGTH; make block length flexible
- Describe installation, configuration, and API
- Get rid of own allocator
