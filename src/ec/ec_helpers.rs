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
                    return P
                }
            }
        }

        pub fn mylogfun<'a>(
            E: &'a Curve,
            P: &'a Point,
            Q: &'a Point,
            k0: usize,
        ) -> impl Fn(&Point) -> (Integer, Integer) + 'a {
            // Compute v = P.weil_pairing(Q, 2^k0)
            // let v = P.weil_pairing(Q, 1 << k0);

            let (v, ok) = E.weil_pairing_2exp(k0 as usize, &P, &Q);
            assert_eq!(ok, 0xFFFFFFFF);

            // Build vs: v, v^2, v^4, ..., v^{2^{k0-1}}
            let mut vs = vec![v.clone()];
            while vs.len() < k0 {
                let last = vs.last().unwrap();
                vs.push(last.square());
            }

            println!("\nmylogfun\n");
            println!("v: {}", v);
            println!();
            // println!("vs: {:?}", vs);
            println!();

            // assert_eq!(vs.last().unwrap(), -&Fq::ONE);
            vs.push(Fq::ONE);
            vs.reverse();

            // Discrete log function
            fn dlog(
                w: &Fq,
                k: usize,
                vs: &[Fq],
            ) -> Integer {
                assert!(k >= 0);
                if k == 0 {
                    // assert!(w.is_one());
                    return 0.big();
                }
                let w2 = w.square();
                let r = dlog(&w2, k - 1, vs);

                let bytes = big_to_bytes(r.clone());
                if vs[k].pow(&bytes, bytes.len() * 8).equals(w) != 0xFFFFFFFF { // TODO: equals
                    // return r + (1 << (k - 1));
                    return r + 2.big().pow(k as u32 - 1)
                }

                r
            }

            // The closure to solve for (a, b)
            move |T: &Point| -> (Integer, Integer) {
                let (wa, ok) = E.weil_pairing_2exp(k0 as usize, &T, &Q);
                assert_eq!(ok, 0xFFFFFFFF);

                println!("wa: {}", wa);
                println!("");

                let (wb, ok) = E.weil_pairing_2exp(k0 as usize, &P, &T);
                assert_eq!(ok, 0xFFFFFFFF);

                println!("wb: {}", wb);
                println!("");

                // let wa = T.weil_pairing(Q, 1 << k0);
                // let wb = P.weil_pairing(T, 1 << k0);
                let a = dlog(&wa, k0, &vs);

                println!("a: {}", a);
                println!("");

                let b = dlog(&wb, k0, &vs);

                println!("b: {}", b);
                println!("");

                // assert!(&(P * Scalar::from(a) + Q * Scalar::from(b)) == T);

                (a, b)
            }
        }
    };
} // End of macro: define_ec_helpers

pub(crate) use define_ec_helpers;
