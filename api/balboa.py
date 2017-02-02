import sys
import os


def _get_env(v):
    _v = os.getenv(v)
    if _v is None:
        sys.stderr.write('Error: Environment variable {0} is undefined\n'.format(v))
        sys.exit(1)
    return _v


def _get_library_suffix():
    if sys.platform == "darwin":
        return 'dylib'
    else:
        return 'so'


def _get_lib_handle(definitions, header, library, build_dir, include_dir):
    from subprocess import Popen, PIPE
    from cffi import FFI
    ffi = FFI()

    interface = Popen(['cc', '-E'] + definitions + [os.path.join(include_dir, header)],
                      stdout=PIPE).communicate()[0].decode("utf-8")
    ffi.cdef(interface)

    suffix = _get_library_suffix()
    lib = ffi.dlopen(os.path.join(build_dir, 'lib', 'lib{0}.{1}'.format(library, suffix)))
    return lib


_build_dir = _get_env('PROJECT_BUILD_DIR')
_include_dir = os.path.dirname(os.path.realpath(__file__))

_lib = _get_lib_handle(
    ['-DBALBOA_API=', '-DBALBOA_NOINCLUDE'],
    'balboa.h',
    'balboa',
    _build_dir,
    _include_dir
)


# outward facing API
new_context = _lib.balboa_new_context
free_context = _lib.balboa_free_context
set_basis = _lib.balboa_set_basis
get_buffer_len = _lib.balboa_get_buffer_len
get_ao = _lib.balboa_get_ao
get_num_aos = _lib.balboa_get_num_aos
