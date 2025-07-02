use criterion::{criterion_group, criterion_main, Criterion};
use klapoti::schemes::klapoti_test::{params64, params128, params256, params512, params768, params1024, params1536, params2048};

pub fn criterion_benchmark(c: &mut Criterion) {
    // c.bench_function("params 64", |b| b.iter(|| params64()));
    // c.bench_function("params 128", |b| b.iter(|| params128()));
    // c.bench_function("params 256", |b| b.iter(|| params256()));
    // c.bench_function("params 512", |b| b.iter(|| params512()));
    // c.bench_function("params 768", |b| b.iter(|| params768()));
    // c.bench_function("params 1024", |b| b.iter(|| params1024()));
    c.bench_function("params 1536", |b| b.iter(|| params1536()));
    // c.bench_function("params 2048", |b| b.iter(|| params2048()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);