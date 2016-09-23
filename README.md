[![Build Status](https://travis-ci.org/bast/balboa.svg?branch=master)](https://travis-ci.org/bast/balboa/builds)

# balboa

Balboa computes Gaussian basis functions and their geometric derivatives for a
batch of points in space.

Experimental code. Under heavy development, nothing is stable.

```
"You know all there is to know about fighting, so there's no sense us going down
that same old road again. To beat this guy, you need speed - you don't have it.
And your knees can't take the pounding, so hard running is out. And you got
arthritis in your neck, and you've got calcium deposits on most of your joints,
so sparring is out.  So, what we'll be calling on is good ol' fashion blunt
force trauma. Horsepower. Heavy-duty, cast-iron, piledriving punches that will
have to hurt so much they'll rattle his ancestors. Every time you hit him with
a shot, it's gotta feel like he tried kissing the express train. Yeah! Let's
start building some hurtin' bombs!" [Rocky Balboa]
```


### Requirements and dependencies

You need the following to install the code:

- C and C++ compilers
- [CMake](https://cmake.org)
- [Python](https://www.python.org)

For testing you need:

- [CFFI](https://cffi.readthedocs.io)
- [Pytest](http://doc.pytest.org)
- [Numpy](http://www.numpy.org)


### Installation

```
virtualenv venv
source venv/bin/activate
pip install -r requirements.txt
./setup
cd build
make
```


### Testing

```
PROJECT_BUILD_DIR=$PWD/build PROJECT_INCLUDE_DIR=$PWD/api PYTHONPATH=$PWD/api py.test -vv test/test.py
```

### TODO

- Get rid of AO_BLOCK_LENGTH; make block length flexible
- Document API


### Ordering of AOs

Let us assume that we have N basis functions and P points.

```
[ geo_000 (undifferentiated)         ][ geo_100 ][ geo_010 ][ geo_001 ] ...
[ ao_1    ][ ao_2    ] ... [ ao_3    ]
[ 1 ... P ][ 1 ... P ] ... [ 1 ... P ]
```
