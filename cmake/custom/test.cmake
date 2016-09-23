enable_testing()

add_test(main_test py.test -vv -s ${PROJECT_SOURCE_DIR}/test/test.py)

set_property(
    TEST
        main_test
    PROPERTY
        ENVIRONMENT PROJECT_BUILD_DIR=${PROJECT_BINARY_DIR}
    )

set_property(
    TEST
        main_test
    APPEND
    PROPERTY
        ENVIRONMENT PROJECT_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/api
    )

set_property(
    TEST
        main_test
    APPEND
    PROPERTY
        ENVIRONMENT PYTHONPATH=${PROJECT_SOURCE_DIR}/api
    )
