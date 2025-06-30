// Fairly ugly-hack benchmark to compute the costs of multiplications, squares
// and inversions for the three base fields. These values are needed to compute
// the correct optimisation costs for the strategies.

use klapoti::fields::Fp512Ext::Fp2 as FpSmall;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::time::Duration;

use klapoti::util::DRNG;

fn m_small(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpSmall = FpSmall::rand(&mut rng);
    let y: FpSmall = FpSmall::rand(&mut rng);

    c.bench_function("New Multiplication (117 bit)", |b| {
        b.iter(|| black_box(x).mul_new(black_box(y)))
    });
}

fn s_small(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpSmall = FpSmall::rand(&mut rng);

    c.bench_function("Square (117 bit)", |b| b.iter(|| black_box(x).square()));
}

fn i_small(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpSmall = FpSmall::rand(&mut rng);

    c.bench_function("Invert (117 bit)", |b| b.iter(|| black_box(x).invert()));
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(3));
    targets = m_small, s_small, i_small
}
criterion_main!(benches);
