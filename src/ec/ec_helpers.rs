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
    };
} // End of macro: define_ec_helpers

pub(crate) use define_ec_helpers;
