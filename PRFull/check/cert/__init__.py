"""PRFull certification library (memory-two friendly rivals, eps -> 0).

    from cert import engine, exact, certify      # with PRFull/check on sys.path

engine  : float and interval leading-order engines (libcert.dylib via ctypes), scans, corner tables
exact   : exact rational engine (gmpy2.mpq or Fraction), E and D*, finite-eps reference
certify : certificate primitives (certify_point, certify_box, certify_ratio_box,
          certify_open_stratum_all_q, exclude_stratum, is_identity)
Build the library with ./compile.sh; run the tests with `nice -n 10 python3 tests.py`.
"""
