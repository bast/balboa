enable_testing()

add_test(test_main py.test -vv -s ${PROJECT_SOURCE_DIR}/test/test.py)

set_property(
    TEST
        test_main
    PROPERTY
        ENVIRONMENT BALBOA_BUILD_DIR=${PROJECT_BINARY_DIR}
    )

set_property(
    TEST
        test_main
    APPEND
    PROPERTY
        ENVIRONMENT PYTHONPATH=${PROJECT_SOURCE_DIR}/balboa
    )

add_test(test_generate py.test -vv -s ${PROJECT_SOURCE_DIR}/balboa/generate.py)
