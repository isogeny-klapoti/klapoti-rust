macro_rules! define_klapoti {
    () => {
        use crate::linalg::matrix::Matrix;
        use crate::quaternion::klpt::klpt;
        use crate::quaternion::lattice::Lattice;
        use crate::quaternion::quadratic_ideal::QuadraticIdeal;
        use crate::quaternion::quadratic_order::{QuadraticOrder, QuadraticOrderEl};
        use crate::quaternion::quaternion_algebra::{QuatAlg, QuatAlgEl};
        use crate::quaternion::quaternion_ideal::QuaternionIdeal;
        use crate::quaternion::quaternion_order::QuaternionOrder;
        use crate::util::{big_to_bytes, bytes_from_str};
        use std::time::Instant;
        use num_traits::Pow;

        /// Let O be an imaginary quadratic order with discriminant D and odd conductor f.
        /// Given an O-oriented supersingular elliptic curve (E, iota), take any omega from O such that O = Z[omega]
        /// and define omega_E := iota(omega).
        /// Let beta from O such that n(omega) + n(beta) = 2^e and gcd(n(beta), n(omega)) = 1. Let P, Q be a basis of E[2^e].
        /// Then the tuple (E, omega, beta, P, Q, omega_E(P), omega_E(Q)) is called a 2dim-representation of (E, iota).
        /// See SCALLOP-HD Definition 9 for more.
        #[derive(Clone, Debug)]
        pub struct TwoDim {
            pub curve: Curve,
            pub omega: QuadraticOrderEl,
            pub beta: QuadraticOrderEl,
            pub P: Point,
            pub Q: Point,
            pub omegaP: Point,
            pub omegaQ: Point,
        }

        fn montgomerize(curve: &Curve, P: &Point) {
            unimplemented!()
        }

        fn canonicalize_orientation(curve: &Curve, P: &Point, Q: &Point, e: u32) -> (Point, Point) {
            let factor = 2.big().pow(e);
            let bytes = big_to_bytes(factor);
            let P4 = curve.mul(&P, &bytes, bytes.len() * 8);
            let Q4 = curve.mul(&Q, &bytes, bytes.len() * 8);

            let PQ4 = curve.add(&P4, &Q4);
            let P2 = curve.mul_small(&P4, 2);
            let Q2 = curve.mul_small(&Q4, 2);

            let points = vec![
                P4.clone(),
                Q4.clone(),
                PQ4.clone(),
                curve.sub(&PQ4, &P2),
                curve.add(&P4, &Q2),
                curve.add(&Q4, &P2),
            ];
            
            for Pt in points.iter() {
                // let iso2 = montgomerize(curve, Pt);
                //let codomain = iso2.codomain();
            }

            // how to use dlog for 2^t:
            let t = 5;
            // let (w1, ok) = curve.weil_pairing_2exp(t, &P, &Q);
            let (w1, ok) = curve.weil_pairing_2exp(246, &P, &Q);
            assert_eq!(ok, 0xFFFFFFFF);


            (*P, *Q)
        }

        impl TwoDim {
            pub fn new(
                curve: Curve,
                omega: QuadraticOrderEl,
                beta: QuadraticOrderEl,
                P: Point,
                Q: Point,
                omegaP: Point,
                omegaQ: Point,
                e: u32,
            ) -> Self {

                canonicalize_orientation(&curve, &P, &Q, e); 

                // replace curve, P, Q,...

                Self {
                    curve,
                    omega,
                    beta,
                    P,
                    Q,
                    omegaP,
                    omegaQ,
                }
            }
        }

        #[derive(Clone, Debug)]
        pub struct PubKey {
            pub product: EllipticProduct,
            pub imagePQ: CouplePoint,
            pub imageOmegaPQ: CouplePoint,
        }

        impl PubKey {
            pub fn new(
                product: EllipticProduct,
                imagePQ: CouplePoint,
                imageOmegaPQ: CouplePoint,
            ) -> Self {
                Self {
                    product,
                    imagePQ,
                    imageOmegaPQ,
                }
            }
        }

        /// KLaPoTi struct
        #[derive(Clone, Debug)]
        pub struct Klapoti {
            pub quadratic_order: QuadraticOrder,
            pub two_dim: TwoDim,
        }

        // TODO: remove
        fn point_order_2e(E: Curve, P: Point, e2: u32) -> u32 {
            let mut T = P.clone();
            for i in 1..=e2 {
                T = E.double(&T);
                if T.isinfinity() == 0xFFFFFFFF {
                    return i; // 2^i
                }
            }
            0 // Not found, not 2^k order
        }

        impl Klapoti {
            pub fn new(quadratic_order: QuadraticOrder, two_dim: TwoDim) -> Self {
                Self {
                    quadratic_order,
                    two_dim,
                }
            }

            pub fn secret(&self) -> QuadraticIdeal {
                self.quadratic_order.random_ideal()
            }

            pub fn act(
                &self,
                ideal: QuadraticIdeal,
                klpt_start_value: u32,
                strategy: Vec<usize>,
            ) -> PubKey {
                let start = Instant::now();

                let disc_abs = self.quadratic_order.order_disc_abs.clone();
                let qa = QuatAlg::new(-disc_abs.clone());

                let e2 = strategy.len() as u32 + 1;

                let basis = Matrix::zeros(4, 4);
                // We use a quadratic order O = Z[(1 + theta)/2].
                // To any pair (beta, gamma) from O^2 we associate a quaternion beta + i * gamma
                // via the identification theta -> j.
                // For this reason we use the quaternion order with basis:
                // [1, i, (1 + j)/2, (i + ij)/2].
                let mut lat = Lattice::new(basis.clone(), 2.big());
                lat.basis[(0, 0)] = 2.big();
                lat.basis[(0, 2)] = 1.big();
                lat.basis[(1, 1)] = 2.big();
                lat.basis[(1, 3)] = 1.big();
                lat.basis[(2, 2)] = 1.big();
                lat.basis[(3, 3)] = 1.big();
                let quaternion_order = QuaternionOrder::new(lat);
                let ideal_norm = ideal.norm();

                let n = ideal.gen1.a.clone();
                assert!(ideal.gen1.b == 0.big());
                let alpha = QuatAlgEl::new(
                    ideal.gen2.a.clone(),
                    0.big(),
                    ideal.gen2.b.clone(),
                    0.big(),
                    ideal.gen2.denom.clone(),
                    qa.clone(),
                );
                let quaternion_ideal = QuaternionIdeal::new_left_ideal(
                    alpha,
                    n.clone(),
                    quaternion_order.clone(),
                    qa.clone(),
                );

                let mut gen_eq = QuatAlgEl::zero(qa.clone());
                let mut found = false;
                loop {
                    for k in klpt_start_value..=e2 {
                        let ok;
                        (ok, gen_eq) = klpt(
                            quaternion_ideal.clone(),
                            qa.clone(),
                            quaternion_order.clone(),
                            k,
                            e2 - 4, // TODO (-3 means as len(strategy))
                        );
                        if ok {
                            found = true;
                            break;
                        }
                    }
                    if found {
                        break;
                    }
                }

                println!("1: {:?}", start.elapsed());
                let second_part = Instant::now();

                gen_eq = gen_eq.normalize();

                let mut gamma_b = QuatAlgEl::new(
                    gen_eq.x.clone(),
                    0.big(),
                    gen_eq.z.clone(),
                    0.big(),
                    gen_eq.denom.clone(),
                    qa.clone(),
                );
                gamma_b = gamma_b.normalize();

                let mut gamma_c = QuatAlgEl::new(
                    gen_eq.y.clone(),
                    0.big(),
                    gen_eq.t.clone(),
                    0.big(),
                    gen_eq.denom.clone(),
                    qa.clone(),
                );
                gamma_c = gamma_c.normalize();

                // TODO: divisions by 2 of gamma_b and gamma_c if needed

                // The two ideals equivalent to the secret ideal `ideal` are then:
                // b = ideal * gamma_b.conj() / norm(ideal)
                // c = ideal * gamma_c.conj() / norm(ideal)

                let norm_b = gamma_b.reduced_norm() / ideal_norm.clone();
                let norm_b = norm_b.numer();

                // The ideal b * c.conj() is principal. The generator is gamma_b.conjugate() * gamma_c / ideal.norm().

                let mut gamma = (gamma_b.conjugate() * gamma_c.clone()) / ideal_norm.clone();
                gamma = gamma.normalize();

                let gamma_quadratic = QuadraticOrderEl::new(
                    gamma.x.clone(),
                    gamma.z.clone(),
                    gamma.denom.clone(),
                    self.quadratic_order.clone(),
                );

                let (u, v) = gamma_quadratic.express_with_el(self.two_dim.omega.clone());

                let u_bytes = big_to_bytes(u);
                let v_bytes = big_to_bytes(v);

                let u_P = self
                    .two_dim
                    .curve
                    .mul(&self.two_dim.P, &u_bytes, u_bytes.len() * 8);
                let u_gammaP =
                    self.two_dim
                        .curve
                        .mul(&self.two_dim.omegaP, &v_bytes, v_bytes.len() * 8);
                let gammaP = self.two_dim.curve.add(&u_P, &u_gammaP);

                let u_Q = self
                    .two_dim
                    .curve
                    .mul(&self.two_dim.Q, &u_bytes, u_bytes.len() * 8);
                let u_gammaQ =
                    self.two_dim
                        .curve
                        .mul(&self.two_dim.omegaQ, &v_bytes, v_bytes.len() * 8);
                let gammaQ = self.two_dim.curve.add(&u_Q, &u_gammaQ);

                let nb_bytes = big_to_bytes(norm_b.clone());

                let norm_b_P =
                    self.two_dim
                        .curve
                        .mul(&self.two_dim.P, &nb_bytes, nb_bytes.len() * 8);
                let norm_b_Q =
                    self.two_dim
                        .curve
                        .mul(&self.two_dim.Q, &nb_bytes, nb_bytes.len() * 8);

                let ell_product = EllipticProduct::new(&self.two_dim.curve, &self.two_dim.curve);

                // TODO
                let fe = 2;
                let PP1 = self.two_dim.curve.mul_small(&norm_b_P, fe);
                let PP2 = self.two_dim.curve.mul_small(&gammaP, fe);

                let QQ1 = self.two_dim.curve.mul_small(&norm_b_Q, fe);
                let QQ2 = self.two_dim.curve.mul_small(&gammaQ, fe);

                
                // debugging:
                let Px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "1226675694967031461374928861727629624321973179136998394931427064315685300027",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "471568703724115901341749099134027127035832042322203601406777353603167167423",
                    )),
                );

                let Py = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "179034482713300109747244485124295104087449802996181620780847970959337093146",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "58543252801583054722328152627939179434526951791836409861241150384998494766",
                    )),
                );

                let PP1 = Point {
                    X: Px,
                    Y: Py,
                    Z: Fq::ONE,
                };


                let Px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "248866086836950554089991015799891182287653190005755700884321151592099782683",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "772889050356731422791785356322718523227084085977542952049919233746522369224",
                    )),
                );

                let Py = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "1161689354733617193117791785489170058103983334851333142158200102813361014218",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "225606592317617417437782069958536757110214474710618758819402131751274917838",
                    )),
                );

                let PP2 = Point {
                    X: Px,
                    Y: Py,
                    Z: Fq::ONE,
                };

                


                let Px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "992421968594890011776270619229553176930598234749523157064938473587440745464",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "969194047278448494096213753293615748421499179500185132821712116215602633182",
                    )),
                );

                let Py = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "1149529860453854205301668921391359438337339035475115004428461471504455462297",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "529273523585589089815020094699351234219495679752062971803541093163272101874",
                    )),
                );

                let QQ1 = Point {
                    X: Px,
                    Y: Py,
                    Z: Fq::ONE,
                };




                let Px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "622853989094175376894443140750699963304765504086282445457642657788184742917",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "474116499916259893879273703802917722069820580712357834431164776605527274420",
                    )),
                );

                let Py = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "450869840280056572265882240833676072037383504959968879400445214195204430674",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "911962211012438493764604134504677168289939912128973267283606877409260575220",
                    )),
                );

                let QQ2 = Point {
                    X: Px,
                    Y: Py,
                    Z: Fq::ONE,
                };

                

                /*
                let P1P2 = CouplePoint::new(&norm_b_P, &gammaP);
                let Q1Q2 = CouplePoint::new(&norm_b_Q, &gammaQ);
                */
                let P1P2 = CouplePoint::new(&PP1, &PP2);
                let Q1Q2 = CouplePoint::new(&QQ1, &QQ2);


                let e22 = 246; // TODO

                let order_foo = point_order_2e(self.two_dim.curve, self.two_dim.P, e22);
                let order_pp1 = point_order_2e(self.two_dim.curve, PP1, e22);
                let order_pp2 = point_order_2e(self.two_dim.curve, PP2, e22);
                let order_qq1 = point_order_2e(self.two_dim.curve, QQ1, e22);
                let order_qq2 = point_order_2e(self.two_dim.curve, QQ2, e22);
  
                println!("");
                println!("");
                println!("order_foo: {}", order_foo);
                println!("order_pp1: {}", order_pp1);
                println!("order_pp2: {}", order_pp2);
                println!("order_qq1: {}", order_qq1);
                println!("order_qq2: {}", order_qq2);
                println!("");




                let inf = Point::INFINITY;

                let image_points = vec![
                    // CouplePoint::new(&self.two_dim.P, &self.two_dim.Q),
                    // CouplePoint::new(&self.two_dim.omegaP, &self.two_dim.omegaQ),
                    CouplePoint::new(&self.two_dim.P, &inf), // (f1(P), f2(P)) or ... ?
                    CouplePoint::new(&self.two_dim.Q, &inf), // (f1(Q), f2(Q)) or ... ?
                    CouplePoint::new(&inf, &self.two_dim.P),
                    CouplePoint::new(&inf, &self.two_dim.Q),
                ];

                println!("");
                println!("-----------");

                println!("");
                println!("norm_b: {:?}", norm_b);

                let norm_c = gamma_c.reduced_norm() / ideal_norm.clone();
                let norm_c = norm_c.numer();
 
                println!("norm_c: {:?}", norm_c);
                println!("");


                let (product, points) = product_isogeny(
                    &ell_product,
                    &P1P2,
                    &Q1Q2,
                    &image_points,
                    e2 as usize,
                    &strategy,
                );

                println!("2: {:?}", second_part.elapsed());

                let goo1 = point_order_2e(product.E1, points[0].P1, e22);

                println!("");
                println!("goo1: {}", goo1);
                println!("");
                println!("");
                println!("");
                println!("");


                let three = Fp::ONE + Fp::ONE + Fp::ONE;
                let four = Fq::ONE + Fq::ONE + Fq::ONE + Fq::ONE;
                let twoh = Fq::new(&Fp::from_i64(256), &Fp::ZERO);

                println!("");
                println!("===============");
                println!("e2: {}", e2);
                println!("");
                println!("E1: {}", product.E1);
                let A = product.E1.A.clone();
                let a2 = A.clone() * A.clone();
                println!("A: {}", a2);
                println!("");
                let num = (a2.clone() - Fq::new(&three, &Fp::ZERO)).clone();
                let num = num.clone() * num.clone() * num; // (A^2 - 3)^3
                let denom = a2 - four;
                let jinv1 = twoh * num / denom;
                println!("");
                println!("E1 j-invariant: {}", jinv1);

                println!("");
                println!("E2: {}", product.E2);
                let A = product.E2.A.clone();
                let a2 = A.clone() * A.clone();
                println!("A: {}", a2);
                println!("");
                let num = (a2.clone() - Fq::new(&three, &Fp::ZERO)).clone();
                let num = num.clone() * num.clone() * num; // (A^2 - 3)^3
                let denom = a2 - four;
                let jinv2 = twoh * num / denom;
                println!("");
                println!("E2 j-invariant: {}", jinv2);

                println!("");
                println!("");

                let ee = 246;
                // TODO: replace hard-coded 246

                let (w1, ok1) = self.two_dim.curve.weil_pairing_2exp(ee, &self.two_dim.P, &self.two_dim.Q);
                assert_eq!(ok1, 0xFFFFFFFF);
                println!("-----------");

                println!("??????????????????????????");
                println!("");
                println!("points[0].P1: {}", &points[0].P1);
                println!("is inf: {}", &points[0].P1.isinfinity());
                println!("");
                println!("points[1].P1: {}", &points[1].P1);
                println!("is inf: {}", &points[1].P1.isinfinity());
                println!("");
                println!("");

                let (w2, ok21) = product.E1.weil_pairing_2exp(ee, &points[0].P1, &points[1].P1);
                // assert_eq!(ok21, 0xFFFFFFFF);
                println!("ok21: {:?}", ok21);

                /*
                let (w2, ok22) = product.E1.weil_pairing_2exp(ee, &points[0].P1, &points[1].P2);
                println!("ok22: {:?}", ok22);

                let (w2, ok23) = product.E1.weil_pairing_2exp(ee, &points[0].P2, &points[1].P1);
                println!("ok23: {:?}", ok23);

                let (w2, ok24) = product.E1.weil_pairing_2exp(ee, &points[0].P2, &points[1].P2);
                println!("ok24: {:?}", ok24);
                */


                let (w22, ok31) = product.E1.weil_pairing_2exp(ee, &points[2].P1, &points[3].P1);
                println!("ok31: {:?}", ok31);

                /*
                let (w22, ok32) = product.E1.weil_pairing_2exp(ee, &points[2].P1, &points[3].P2);
                println!("ok32: {:?}", ok32);

                let (w22, ok33) = product.E1.weil_pairing_2exp(ee, &points[2].P2, &points[3].P1);
                println!("ok33: {:?}", ok33);

                let (w22, ok34) = product.E1.weil_pairing_2exp(ee, &points[2].P2, &points[3].P2);
                println!("ok34: {:?}", ok34);
                */

                // let norm_b_u32 = norm_b.to_u32_wrapping();
                let bytes1 = big_to_bytes(norm_b.clone());
                let foo1 = w2.pow(&bytes1, bytes1.len() * 8);
 
                let bytes2 = big_to_bytes(norm_c.clone());
                let foo2 = w2.pow(&bytes2, bytes2.len() * 8);

                println!("");
                println!("bytes1: {:?}", bytes1);
                println!("");

                println!("w1: {:?}", w1);
                println!("w1: {}", w1);
                println!("");

                // Why is w2 = 1 ????????????????
                println!("w2: {:?}", w2);
                println!("w2: {}", w2);
                println!("");

                println!("w22: {:?}", w22);
                println!("w22: {}", w22);
                println!("");

                println!("foo1: {:?}", foo1);
                println!("");
                println!("foo2: {:?}", foo2);
                println!("");
                println!("--------");


                PubKey::new(product, points[0], points[1])
            }
        }
    };
} // End of macro: define_klapoti

pub(crate) use define_klapoti;
