#[macro_use]
extern crate criterion;
extern crate least_action;

use criterion::Criterion;
use criterion::black_box;
use least_action::*;
use std::iter;

fn free_body_1d_brute_force(n: i64, m: usize) -> Vec<i64> {
    let p = Phys1D::new(
        vec![],
        free_body,
        0,
        n+1,
        1f64
    );

    p.brute_force(m)
}

fn bench_fb_1d_bf(c: &mut Criterion) {
    c.bench_function_over_inputs(
        "Freebody 1d bruteforce",
        |b, size: &i64 | {
            b.iter(|| free_body_1d_brute_force(*size, 2))
        },
        vec![50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
    );
}

criterion_group!(benches, bench_fb_1d_bf);
criterion_main!(benches);