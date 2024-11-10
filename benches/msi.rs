// Fairly ugly-hack benchmark to compute the costs of multiplications, squares
// and inversions for the three base fields. These values are needed to compute
// the correct optimisation costs for the strategies.

use klapoti::fields::Fp117Ext::Fp2 as FpSmall;
use klapoti::fields::Fp1757Ext::Fp2 as FpBig;
use klapoti::fields::Fp214Ext::Fp2 as FpMid;
use klapoti::fields::Fp509Ext::Fp2 as FpMiddle;

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

fn m_mid(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpMid = FpMid::rand(&mut rng);
    let y: FpMid = FpMid::rand(&mut rng);

    c.bench_function("New Multiplication (214 bit)", |b| {
        b.iter(|| black_box(x).mul_new(black_box(y)))
    });
}

fn s_mid(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpMid = FpMid::rand(&mut rng);

    c.bench_function("Square (214 bit)", |b| b.iter(|| black_box(x).square()));
}

fn i_mid(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpMid = FpMid::rand(&mut rng);

    c.bench_function("Invert (214 bit)", |b| b.iter(|| black_box(x).invert()));
}

fn m_middle(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpMiddle = FpMiddle::rand(&mut rng);
    let y: FpMiddle = FpMiddle::rand(&mut rng);

    c.bench_function("New Multiplication (509 bit)", |b| {
        b.iter(|| black_box(x).mul_new(black_box(y)))
    });
}

fn s_middle(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpMiddle = FpMiddle::rand(&mut rng);

    c.bench_function("Square (509 bit)", |b| b.iter(|| black_box(x).square()));
}

fn i_middle(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpMiddle = FpMiddle::rand(&mut rng);

    c.bench_function("Invert (509 bit)", |b| b.iter(|| black_box(x).invert()));
}

fn m_big(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpBig = FpBig::rand(&mut rng);
    let y: FpBig = FpBig::rand(&mut rng);

    c.bench_function("New Multiplication (1745 bit)", |b| {
        b.iter(|| black_box(x).mul_new(black_box(y)))
    });
}

fn s_big(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpBig = FpBig::rand(&mut rng);

    c.bench_function("Square (1745 bit)", |b| b.iter(|| black_box(x).square()));
}

fn i_big(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: FpBig = FpBig::rand(&mut rng);

    c.bench_function("Invert (1745 bit)", |b| b.iter(|| black_box(x).invert()));
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(3));
    targets = m_small, s_small, i_small, m_big, s_big, i_big, m_mid, s_mid, i_mid, m_middle, s_middle, i_middle
}
criterion_main!(benches);
