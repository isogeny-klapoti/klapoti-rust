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
        use crate::util::{big_to_bytes, valuation, bytes_from_str};
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
                P: Point,
                Q: Point,
                omegaP: Point,
                omegaQ: Point,
                valuation_2: u32, // TODO: move
                cofactor: u32, // TODO: move
            ) -> Self {

                // TODO:
                curve.normalize();

                let bytes = big_to_bytes(2.big().pow(valuation_2 - 1));

                let R = generate_random_fq(&curve, (valuation_2 - 1).big(), cofactor.big());
                let mut S = Point::INFINITY; 
                loop {
                    S = generate_random_fq(&curve, (valuation_2 - 1).big(), cofactor.big());
                    let (w, ok) = curve.weil_pairing_2exp(valuation_2 as usize, &R, &S);
                    assert_eq!(ok, 0xFFFFFFFF);
                    let wto = w.pow(&bytes, bytes.len() * 8);

                    if wto.equals(&Fq::ONE) == 0 { // wto != 1
                        break;
                    } 
                }


                let A = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "87588114400902199998338770337694739636315554156424733587465955178621864700",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "619983003027140970448668231864152961217819436005861629692568717681019587071",
                    )),
                );
                let curve = Curve::new(&A);

                let Px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "640262070225397919641419391798411906185192704400982985445972594570331098653",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "287810177297146274761114922305202226260110532495201421327960422009294927864",
                    )),
                );
                let Py = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "457788924259115902463588305888251116114261190788876950068401387037922249228",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "633138429939111716122288253451612328751828562769273084337837985144875427010",
                    )),
                );
                let P = Point {
                    X: Px,
                    Y: Py,
                    Z: Fq::ONE,
                };

                let Px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "577933924249040245437036418717899949221523555485972514647508407144730672525",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "1019387042211059007922076628149515744324272062266127407040430254585022649375",
                    )),
                );
                let Py = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "509394661907087766331273068402842237724830983103433527961227534947863445264",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "725385452187303835136549965145606007505912062894778454060531039947591434989",
                    )),
                );
                let Q = Point {
                    X: Px,
                    Y: Py,
                    Z: Fq::ONE,
                };


                // omegaP, omegaQ

                let Px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "634039621850969308805385027694276169360651609855271563211520020735844785215",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "1075188422265331762811290023078590423058046429393430765873987622134783102482",
                    )),
                );
                let Py = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "1097972323119634683247217098336541543090612946284328422947693398582605586101",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "220599266309031818624969048679582813355156390507263420502962613496935193172",
                    )),
                );
                let omegaP = Point {
                    X: Px,
                    Y: Py,
                    Z: Fq::ONE,
                };

                let Px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "986992327121915758581192528730652032676947074727662307563723961364053414999",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "921059537265102690615221739730873763647293364847969554649992493452361646901",
                    )),
                );
                let Py = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "74607066523784015036897537654048456971960892140972951534399235879024055198",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "371448068294077661734717191164184817605554068703830199822253857365493561338",
                    )),
                );
                let omegaQ = Point {
                    X: Px,
                    Y: Py,
                    Z: Fq::ONE,
                };

                /*
                println!("");
                println!("+++++++++++++++++++++++++++++++++");
                println!("");
                println!("R: {}", R);
                println!("");
                println!("S: {}", S);
                println!("");
                println!("");
                */

                let Px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "925473869824409399843652150961237538577304088344928976475296163379511372554",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "116287816262904234454738007257879993288471265874250810570974463204167683390",
                    )),
                );
                let Py = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "1015993309849213971252813187901840334758687159880946344187583207756725238062",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "957009723818330665698232654864907677082438236295236482958661892033234280556",
                    )),
                );
                let R = Point {
                    X: Px,
                    Y: Py,
                    Z: Fq::ONE,
                };

                let Px = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "838916474589319857260876467036924575254055493201918407630183323525256184301",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "1116028115080817570187543946914933769013649727010082549357974164345300672012",
                    )),
                );
                let Py = Fq::new(
                    &Fp::decode_reduce(&bytes_from_str(
                        "227413930184305938291645565787496783486997670939025411806399471552837775354",
                    )),
                    &Fp::decode_reduce(&bytes_from_str(
                        "69513907430609967024884397322368381340309250746433465152060923612821002173",
                    )),
                );
                let S = Point {
                    X: Px,
                    Y: Py,
                    Z: Fq::ONE,
                };

                let mylog = mylogfun(&curve, &P, &Q, valuation_2 as usize);

                println!("");
                println!("--------------------");
                println!("");

                let roo = mylog(&R);
                let soo = mylog(&S);

                let a = roo.0.clone();
                let b = roo.1.clone();
                let c = soo.0.clone();
                let d = soo.1.clone();

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

                

                println!("");
                println!("m00: {}", m00);
                println!("");
                println!("m01: {}", m01);
                println!("");
                println!("m10: {}", m10);
                println!("");
                println!("m11: {}", m11);
                println!("");
                println!("");
                println!("");

                /*
                let m00 = ( d * inv_det).rem_euclid(p);
                let m01 = ((-b) * inv_det).rem_euclid(p);
                let m10 = ((-c) * inv_det).rem_euclid(p);
                let m11 = ( a * inv_det).rem_euclid(p);
                */

                let mut mat_inv = Matrix::<Integer>::zeros(4, 4);
                mat_inv[(0, 0)] = m00;
                mat_inv[(0, 1)] = m01;
                mat_inv[(1, 0)] = m10;
                mat_inv[(1, 1)] = m11;

                let foo = mat * mat_inv;
                println!("");
                println!("foo: {}", foo[(0, 0)].clone().modulo(&m));
                println!("");
                println!("foo: {}", foo[(0, 1)].clone().modulo(&m));
                println!("");
                println!("foo: {}", foo[(1, 0)].clone().modulo(&m));
                println!("");
                println!("foo: {}", foo[(1, 1)].clone().modulo(&m));
                println!("");

                let omega_roo = mylog(&omegaP);
                let omega_soo = mylog(&omegaQ);

                println!("");
                println!("roo: {}, {}", roo.0, roo.1);
                println!("");
                println!("soo: {}, {}", soo.0, soo.1);
                println!("");
                println!("");

                println!("omega_roo: {}, {}", omega_roo.0, omega_roo.1);
                println!("");
                println!("omega_soo: {}, {}", omega_soo.0, omega_soo.1);
                println!("");
                println!("");


                Self {
                    curve,
                    P,
                    Q,
                    omegaP,
                    omegaQ,
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

                let norm_b = "4565924146296632737235756772431642587386221405171208750233525117839286963".big();
                let norm_c = "2501464112816904581097433230540031475923714182331267082252899687331192141".big();

                // The ideal b * c.conj() is principal. The generator is gamma_b.conjugate() * gamma_c / ideal.norm().

                let mut gamma = (gamma_b.conjugate() * gamma_c.clone()) / ideal_norm.clone();
                gamma = gamma.normalize();

                // debugging:
                // γ = -298039931075527522255844797466619404933784424563074793597744939*ϑ - 3127764656277006504647309206099010060132678986748243084380759488225833242

                let mut gamma = QuatAlgEl::new(
                    "-3127764656277006504647309206099010060132678986748243084380759488225833242".big(),
                    0.big(),
                    "-298039931075527522255844797466619404933784424563074793597744939".big(),
                    0.big(),
                    1.big(),
                    qa.clone(),
                );


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
                println!("e_start: {}", e_start);
                println!("");


                let fe = valuation_2 - 2 - e_start;
                let fe = 2.big().pow(fe);
                let fe_bytes = big_to_bytes(fe);

                let fe = 4; // TODO: hardcoded for debugging
                let fe_bytes = big_to_bytes(fe.big()); // TODO: hardcoded for debugging

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

                println!("");
                println!("key: {}, {}", e, e-1);
                println!("");

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

                println!("");
                println!("e: {}", e);
                println!("");

                let three = Fp::ONE + Fp::ONE + Fp::ONE;
                let four = Fq::ONE + Fq::ONE + Fq::ONE + Fq::ONE;
                let twoh = Fq::new(&Fp::from_i64(256), &Fp::ZERO);

                println!("");
                println!("===============");
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

                let (z, ok1) = self.two_dim.curve.weil_pairing_2exp(valuation_2 as usize, &self.two_dim.P, &self.two_dim.Q);
                assert_eq!(ok1, 0xFFFFFFFF);

                let ztob = z.pow(&nb_bytes, nb_bytes.len() * 8);

                let distinguish = || -> (usize, Point, Point) {
                    let imP1 = vec![&points[0].P1, &points[0].P2];
                    let imQ1 = vec![&points[1].P1, &points[1].P2];

                    // let imP2 = vec![&points[2].P1, &points[2].P2]; // TODO: the corresponding points were removed
                    // let imQ2 = vec![&points[3].P1, &points[3].P2];

                    println!("");
                    println!("============ ev ==============");
                    println!("");
                    
                    for i in 0..2 {
                        let imP = &imP1[i];
                        let imQ = &imQ1[i];
                        let mut curve = product.E1;
                        if i == 1 {
                            curve = product.E2;
                        }

                        println!("--------- 11");
                        println!("imP1: {}, {}", imP.X / imP.Z, imP.Y / imP.Z);
                        println!("");
                        println!("imQ1: {}, {}", imQ.X / imQ.Z, imQ.Y / imQ.Z);
                        println!("");
                        println!("");

                        let (w, ok) = curve.weil_pairing_2exp(valuation_2 as usize, &imP, &imQ);
                        println!("ok: {:?}", ok);
                        println!("");
                        println!("z^Nb: {}", ztob);
                        println!("");
                        println!("w: {}", w);
                        println!("");
                        let w_inv = w.invert();
                        println!("w inv: {}", w_inv);
                        println!("");

                        if ztob == w || ztob == w_inv {
                            println!("Found match for ztob!");
                            if i == 0 {
                                // norm_b "is" curve 1
                            } else {
                                // norm_b "is" curve 2
                            }
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
                println!("ok 1: {:?}", ok);
                if ok == 0 || w != ztob {
                    imQ.set_neg();
                }

                let norm_omega = self.two_dim.omega.norm();
                println!("norm_omega: {}", norm_omega);

                let bytes1 = big_to_bytes(norm_b * norm_omega.clone());
                let ztow = z.pow(&bytes1, bytes1.len() * 8);
                println!("ztow: {}", ztow);

                let im_omegaP12 = vec![&points[2].P1, &points[2].P2];
                let im_omegaQ12 = vec![&points[3].P1, &points[3].P2];

                let im_omegaP = im_omegaP12[ind];
                let im_omegaQ = im_omegaQ12[ind];

                let (c2, ok) = curve.weil_pairing_2exp(valuation_2 as usize, &im_omegaP, &im_omegaQ);
                println!("ok 3: {:?}", ok);
                if ok == 0 || c2 != ztow {
                    imQ.set_neg();
                }
                

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

                PubKey::new(curve, imP, imQ, *im_omegaP, *im_omegaQ, valuation_2, cofactor)
            }
        }
    };
} // End of macro: define_klapoti

pub(crate) use define_klapoti;
