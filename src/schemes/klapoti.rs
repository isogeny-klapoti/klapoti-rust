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
        use crate::util::{big_to_bytes, valuation};
        use std::time::Instant;
        use num_traits::Pow;
        use std::collections::HashMap;

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

        impl TwoDim {
            pub fn new(
                curve: Curve,
                omega: QuadraticOrderEl,
                beta: QuadraticOrderEl,
                P: Point,
                Q: Point,
                omegaP: Point,
                omegaQ: Point,
            ) -> Self {
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
            pub curve: Curve,
            pub P: Point,
            pub Q: Point,
            pub omegaP: Point,
            pub omegaQ: Point,
        }

        impl PubKey {
            pub fn new(
                curve: Curve,
                mut P: Point,
                mut Q: Point,
                mut omegaP: Point,
                mut omegaQ: Point,
                valuation_2: u32,
                cofactor: u32,
            ) -> Self {

                let (new_curve, isom) = curve.normalize();
                new_curve.ec_iso_eval(&mut P, &isom);
                new_curve.ec_iso_eval(&mut Q, &isom);
                new_curve.ec_iso_eval(&mut omegaP, &isom);
                new_curve.ec_iso_eval(&mut omegaQ, &isom);

                let Px = PointX::new_xz(&P.X, &P.Z);
                let (P, _) = new_curve.complete_pointX(&Px);

                let Qx = PointX::new_xz(&Q.X, &Q.Z);
                let (Q, _) = new_curve.complete_pointX(&Qx);

                let omegaPx = PointX::new_xz(&omegaP.X, &omegaP.Z);
                let (omegaP, _) = new_curve.complete_pointX(&omegaPx);

                let omegaQx = PointX::new_xz(&omegaQ.X, &omegaQ.Z);
                let (omegaQ, _) = new_curve.complete_pointX(&omegaQx);

                let bytes = big_to_bytes(2.big().pow(valuation_2 - 1));
                let R = generate_random_fq(&new_curve, (valuation_2 - 1).big(), cofactor.big());
                
                let mut S;
                loop {
                    S = generate_random_fq(&new_curve, (valuation_2 - 1).big(), cofactor.big());
                    let (w, ok) = new_curve.weil_pairing_2exp(valuation_2 as usize, &R, &S);
                    assert_eq!(ok, 0xFFFFFFFF);
                    let wto = w.pow(&bytes, bytes.len() * 8);

                    if wto.equals(&Fq::ONE) == 0 { // wto != 1
                        break;
                    } 
                }
                
                let dlog = prepare_dlog_solver(&new_curve, &P, &Q, valuation_2 as usize);
                let rdlog = dlog(&R);
                let sdlog = dlog(&S);

                let a = rdlog.0.clone();
                let b = rdlog.1.clone();
                let c = sdlog.0.clone();
                let d = sdlog.1.clone();

                let mut mat = Matrix::<Integer>::zeros(4, 4);
                mat[(0, 0)] = a.clone();
                mat[(0, 1)] = b.clone();
                mat[(1, 0)] = c.clone();
                mat[(1, 1)] = d.clone();

                let m =  2.big().pow(valuation_2);
                let det = (a.clone() * d.clone() - b.clone() * c.clone()).modulo(&m);
                let det_inv = det.invert(&m).unwrap();
                let m00 = (d * det_inv.clone()).modulo(&m);
                let m01 = ((-b).modulo(&m) * det_inv.clone()).modulo(&m);
                let m10 = ((-c).modulo(&m) * det_inv.clone()).modulo(&m);
                let m11 = (a * det_inv).modulo(&m);

                let mut mat_inv = Matrix::<Integer>::zeros(4, 4);
                mat_inv[(0, 0)] = m00;
                mat_inv[(0, 1)] = m01;
                mat_inv[(1, 0)] = m10;
                mat_inv[(1, 1)] = m11;
                
                let omega_rdlog = dlog(&omegaP);
                let omega_dlog = dlog(&omegaQ);

                let mut mat_om = Matrix::<Integer>::zeros(4, 4);
                mat_om[(0, 0)] = omega_rdlog.0.clone();
                mat_om[(0, 1)] = omega_rdlog.1.clone();
                mat_om[(1, 0)] = omega_dlog.0.clone();
                mat_om[(1, 1)] = omega_dlog.1.clone();

                let omega_RS = mat * mat_om * mat_inv;

                let mut bytes = big_to_bytes(omega_RS[(0, 0)].clone());
                let mut T1 = new_curve.mul(&R, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(omega_RS[(0, 1)].clone());
                let mut T2 = new_curve.mul(&S, &bytes, bytes.len() * 8);
                let omegaR = new_curve.add(&T1, &T2);

                let mut bytes = big_to_bytes(omega_RS[(1, 0)].clone());
                T1 = new_curve.mul(&R, &bytes, bytes.len() * 8);
                bytes = big_to_bytes(omega_RS[(1, 1)].clone());
                T2 = new_curve.mul(&S, &bytes, bytes.len() * 8);
                let omegaS = new_curve.add(&T1, &T2);

                Self {
                    curve: new_curve,
                    P: R,
                    Q: S,
                    omegaP: omegaR,
                    omegaQ: omegaS,
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
                strategies: HashMap<u32, Vec<usize>>,
                valuation_2: u32,
                cofactor: u32, // p + 1 = 2^valuation_2 * cofactor
            ) -> PubKey {
                let start = Instant::now();

                let disc_abs = self.quadratic_order.order_disc_abs.clone();
                let qa = QuatAlg::new(-disc_abs.clone());

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
                    for k in klpt_start_value..=valuation_2 {
                        let ok;
                        (ok, gen_eq) = klpt(
                            quaternion_ideal.clone(),
                            qa.clone(),
                            quaternion_order.clone(),
                            k,
                            valuation_2 - 4, // - 2 - 2 because there is 2^2 factor in KLPT
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

                // The two ideals equivalent to the secret ideal `ideal` are then:
                // b = ideal * gamma_b.conj() / norm(ideal)
                // c = ideal * gamma_c.conj() / norm(ideal)

                let norm_b = gamma_b.reduced_norm() / ideal_norm.clone();
                let norm_b = norm_b.numer();

                let norm_c = gamma_c.reduced_norm() / ideal_norm.clone();
                let norm_c = norm_c.numer();

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
                let u_bytes = big_to_bytes(u.clone());
                let v_bytes = big_to_bytes(v.clone());

                let mut u_P = self
                    .two_dim
                    .curve
                    .mul(&self.two_dim.P, &u_bytes, u_bytes.len() * 8);
                if u < 0.big() {
                    u_P.set_neg();
                }

                let mut v_omegaP =
                    self.two_dim
                        .curve
                        .mul(&self.two_dim.omegaP, &v_bytes, v_bytes.len() * 8);
                if v < 0.big() {
                    v_omegaP.set_neg();
                }
                let gammaP = self.two_dim.curve.add(&u_P, &v_omegaP);

                let mut u_Q = self
                    .two_dim
                    .curve
                    .mul(&self.two_dim.Q, &u_bytes, u_bytes.len() * 8);
                if u < 0.big() {
                    u_Q.set_neg();
                }
                let mut v_omegaQ =
                    self.two_dim
                        .curve
                        .mul(&self.two_dim.omegaQ, &v_bytes, v_bytes.len() * 8);
                if v < 0.big() {
                    v_omegaQ.set_neg();
                }
                let gammaQ = self.two_dim.curve.add(&u_Q, &v_omegaQ);

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

                let e_start = valuation(Integer::from(norm_b.clone()) + Integer::from(norm_c.clone()), Integer::from(2)).0 as u32;

                let fe = valuation_2 - 2 - e_start;
                let fe = 2.big().pow(fe);
                let fe_bytes = big_to_bytes(fe);

                let mut PP1 = self.two_dim.curve.mul(&norm_b_P, &fe_bytes, fe_bytes.len() * 8);
                let mut PP2 = self.two_dim.curve.mul(&gammaP, &fe_bytes, fe_bytes.len() * 8);

                let mut QQ1 = self.two_dim.curve.mul(&norm_b_Q, &fe_bytes, fe_bytes.len() * 8);
                let mut QQ2 = self.two_dim.curve.mul(&gammaQ, &fe_bytes, fe_bytes.len() * 8);

                fn h2(curve: &Curve, T1: &Point, T2: &Point) -> (Point, Point) {
                    (curve.add(&T1, &T2), curve.sub(&T1, &T2))
                }

                let mut e = e_start;
                loop {
                    let mut T = self.two_dim.curve.sub(&PP1, &PP2);
                    for _ in 0..e+1 {
                        T = self.two_dim.curve.double(&T);
                    }
                    // T is now 2^(e+1) * (PP1 - PP2)
                    if T.isinfinity() == 0xFFFFFFFF {
                        (PP1, PP2) = h2(&self.two_dim.curve, &PP1, &PP2);
                        (QQ1, QQ2) = h2(&self.two_dim.curve, &QQ1, &QQ2);
                        e -= 1;
                    } else {
                        break;
                    }
                }

                let P1P2 = CouplePoint::new(&PP1, &PP2);
                let Q1Q2 = CouplePoint::new(&QQ1, &QQ2);

                let order_foo = point_order_2e(self.two_dim.curve, self.two_dim.P, valuation_2);
                let order_pp1 = point_order_2e(self.two_dim.curve, PP1, valuation_2);

                println!("");
                println!("");
                println!("P: {}, {}", self.two_dim.P.X / self.two_dim.P.Z, self.two_dim.P.Y / self.two_dim.P.Z);
                println!("");
                println!("valuation_2: {}", valuation_2);
                println!("");
                println!("curve: {}", self.two_dim.curve);
                println!("");
                println!("order_foo: {}", order_foo);
                println!("order_pp1: {}", order_pp1);
                println!("");

                let pre = |T1: &Point, T2: &Point| -> (Point, Point){
                    let mut K1 = T1.clone();
                    let mut K2 = T2.clone();
                    for _ in 0..(e_start - e) {
                        let (new_T1, new_T2) = h2(&self.two_dim.curve, &K1, &K2);
                        K1 = new_T1;
                        K2 = new_T2;
                    }
                
                    (K1, K2)
                };

                let inf = Point::INFINITY;
                let image_points = vec![
                    {
                        let (p1, p2) = pre(&self.two_dim.P, &inf);
                        CouplePoint::new(&p1, &p2)
                    },
                    {
                        let (p1, p2) = pre(&self.two_dim.Q, &inf);
                        CouplePoint::new(&p1, &p2)
                    },
                    {
                        let (p1, p2) = pre(&self.two_dim.omegaP, &inf);
                        CouplePoint::new(&p1, &p2)
                    },
                    {
                        let (p1, p2) = pre(&self.two_dim.omegaQ, &inf);
                        CouplePoint::new(&p1, &p2)
                    },
                ];

                if strategies.get(&(e-1)).is_none() {
                    panic!("No strategy for e - 1 = {}", e - 1);
                }
                
                let (product, points) = product_isogeny(
                    &ell_product,
                    &P1P2,
                    &Q1Q2,
                    &image_points,
                    e as usize,
                    &strategies[&(e-1)],
                );

                println!("2: {:?}", second_part.elapsed());
                let third_part = Instant::now();

                let distinguish_part = Instant::now();

                let (z, ok1) = self.two_dim.curve.weil_pairing_2exp(valuation_2 as usize, &self.two_dim.P, &self.two_dim.Q);
                assert_eq!(ok1, 0xFFFFFFFF);

                let ztob = z.pow(&nb_bytes, nb_bytes.len() * 8);

                let distinguish = || -> (usize, Point, Point) {
                    let imP1 = vec![&points[0].P1, &points[0].P2];
                    let imQ1 = vec![&points[1].P1, &points[1].P2];

                    for i in 0..2 {
                        let imP = &imP1[i];
                        let imQ = &imQ1[i];
                        let mut curve = product.E1;
                        if i == 1 {
                            curve = product.E2;
                        }

                        let (w, ok) = curve.weil_pairing_2exp(valuation_2 as usize, &imP, &imQ);
                        assert_eq!(ok, 0xFFFFFFFF);
                        let w_inv = w.invert();

                        if ztob == w || ztob == w_inv {
                            return (i, *imP.clone(), *imQ.clone());
                        }
                    }
                    return (2, Point::INFINITY, Point::INFINITY);
                };

                let (ind, imP, mut imQ) = distinguish();
                assert!(ind != 2);
                let mut curve = product.E1;
                if ind == 1 {
                    curve = product.E2;
                }

                let (w, ok) = curve.weil_pairing_2exp(valuation_2 as usize, &imP, &imQ);
                if ok == 0 || w != ztob {
                    imQ.set_neg();
                }

                let norm_omega = self.two_dim.omega.norm();
                let bytes1 = big_to_bytes(norm_b * norm_omega.clone());
                let ztow = z.pow(&bytes1, bytes1.len() * 8);

                let im_omegaP12 = vec![&points[2].P1, &points[2].P2];
                let im_omegaQ12 = vec![&points[3].P1, &points[3].P2];
                let im_omegaP = im_omegaP12[ind];
                let im_omegaQ = im_omegaQ12[ind];

                let (c2, ok) = curve.weil_pairing_2exp(valuation_2 as usize, &im_omegaP, &im_omegaQ);
                if ok == 0 || c2 != ztow {
                    imQ.set_neg();
                }

                println!("distinguish: {:?}", distinguish_part.elapsed());
                

                println!("");
                println!("???????????????? ============");
                println!("");
                println!("imP: {}, {}", imP.X / imP.Z, imP.Y / imP.Z);
                println!("");
                println!("imQ: {}, {}", imQ.X / imQ.Z, imQ.Y / imQ.Z);
                println!("");

                println!("im_omegaP: {}, {}", im_omegaP.X / im_omegaP.Z, im_omegaP.Y / im_omegaP.Z);
                println!("");
                println!("im_omegaQ: {}, {}", im_omegaQ.X / im_omegaQ.Z, im_omegaQ.Y / im_omegaQ.Z);
                println!("");
                println!("curve: {}", curve);
                println!("");

                let pub_key = PubKey::new(curve, imP, imQ, *im_omegaP, *im_omegaQ, valuation_2, cofactor);

                println!("");
                println!("");
                println!("");
                println!("3: {:?}", third_part.elapsed());
                println!("");

                pub_key
            }
        }
    };
} // End of macro: define_klapoti

pub(crate) use define_klapoti;
