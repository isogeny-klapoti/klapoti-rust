#![allow(non_snake_case)]

macro_rules! define_ec_helpers {
    () => {
        use rand::prelude::*;
        use rand_chacha::ChaCha20Rng;

        /// Generate a random point on the curve with x in Fq.
        /// The point needs to be of order `x^k`, where `order_minus = x^{k-1}`.
        pub fn generate_random_fq(curve: &Curve, order_minus: Integer, f: Integer) -> Point {
            // p + 1 = order * f
            let mut rng = ChaCha20Rng::from_entropy();
            let mut X: Fq;
            let mut P: Point;
            let f_bytes = big_to_bytes(f);

            loop {
                X = Fq::rand(&mut rng);
                let Pxz = PointX::new_xz(&X, &Fq::ONE);
                let (Px, _) = curve.complete_pointX(&Pxz);
                P = curve.mul(&Px, &f_bytes, f_bytes.len() * 8);

                let bytes = big_to_bytes(order_minus.clone());
                let Q = curve.mul(&P, &bytes, bytes.len() * 8);

                if Q.isinfinity() == 0x00000000 {
                    return P;
                }
            }
        }

        pub fn prepare_dlog_solver<'a>(
            E: &'a Curve,
            P: &'a Point,
            Q: &'a Point,
            k0: usize,
        ) -> impl Fn(&Point) -> (Integer, Integer) + 'a {
            let (v, ok) = E.weil_pairing_2exp(k0 as usize, &P, &Q);
            assert_eq!(ok, 0xFFFFFFFF);

            // Build vs: v, v^2, v^4, ..., v^{2^{k0-1}}
            let mut vs = vec![v.clone()];
            while vs.len() < k0 {
                let last = vs.last().unwrap();
                vs.push(last.square());
            }

            // Discrete log function
            fn dlog(w: &Fq, k: usize, vs: &[Fq], k0: usize) -> Integer {
                if k == 0 {
                    return 0.big();
                }
                if k == 1 {
                    return if w.equals(&Fq::ONE) == 0xFFFFFFFF { 0.big() } else { 1.big() };
                }

                let m = (k + 1) >> 1;
                let l = k - m;

                let mut wm = w.clone();
                for _ in 0..m {
                    wm.set_square();
                }
                let y = dlog(&wm, l, vs, k0);

                let ybytes = big_to_bytes(y.clone());
                let mut vy = vs[k0-k].clone();
                vy.set_pow_simple(&ybytes);
                let wl = w / vy;
                let z = dlog(&wl, m, vs, k0);

                return y + 2.big().pow(l as u32) * z;
            }

            // The closure to solve for (a, b)
            move |T: &Point| -> (Integer, Integer) {
                let (wa, ok) = E.weil_pairing_2exp(k0 as usize, &T, &Q);
                assert_eq!(ok, 0xFFFFFFFF);

                let (wb, ok) = E.weil_pairing_2exp(k0 as usize, &P, &T);
                assert_eq!(ok, 0xFFFFFFFF);

                let a = dlog(&wa, k0, &vs, k0);
                let b = dlog(&wb, k0, &vs, k0);

//                let abytes = big_to_bytes(a.clone());
//                let bbytes = big_to_bytes(b.clone());
//                assert!(E.add(&E.mul(&P, &abytes, 8*abytes.len()), &E.mul(&Q, &bbytes, 8*bbytes.len())).equals(T) == 0xFFFFFFFF);

                (a, b)
            }
        }
    };
} // End of macro: define_ec_helpers

pub(crate) use define_ec_helpers;
