# README

This work uses the library
[two-isogenies](https://github.com/ThetaIsogenies/two-isogenies) for the computation of (2, 2)-isogenies.
To use its function `product_isogeny` a strategy has to be precomputed using the script `strategy.py`:

```
data = [1733, 2937, 2355, 85162]
strat = optimised_strategy(*data)
print(strat)
```

Note that the list `data` contains the isogeny chain, the cost of multiplication, the cost of squaring,
and the cost of inversion. These costs are computed using a benchmark `benches/msi.rs`, just run:
```
cargo bench
```

Example output:
```
New Multiplication (1745 bit)
                        time:   [2.9330 µs 2.9375 µs 2.9424 µs]
                        change: [-2.4460% -1.9078% -1.4550%] (p = 0.00 < 0.05)
                        Performance has improved.
Found 3 outliers among 100 measurements (3.00%)
  2 (2.00%) high mild
  1 (1.00%) high severe

Square (1745 bit)       time:   [2.3517 µs 2.3552 µs 2.3589 µs]
                        change: [-1.8670% -1.6079% -1.3693%] (p = 0.00 < 0.05)
                        Performance has improved.
Found 13 outliers among 100 measurements (13.00%)
  6 (6.00%) low mild
  4 (4.00%) high mild
  3 (3.00%) high severe

Invert (1745 bit)       time:   [84.972 µs 85.162 µs 85.382 µs]
                        change: [-2.8793% -2.6136% -2.3365%] (p = 0.00 < 0.05)

```

As costs we take the middle values of the `time` list.

## Tests

To run the tests:
```
cargo test
```

To run a test, execute for example:
```
cargo test params512
```

To run a test in a release mode (much faster):
```
cargo test --release params512
```

Or with print statements:
```
cargo test --release -- --nocapture params512
```