#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

pub mod ec;
mod error;
pub mod fields;
pub mod finitefield;
pub mod linalg;
pub mod quaternion;
pub mod schemes;
mod theta;
pub mod util;

pub mod ec117 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp117::Fp;
    pub type Fq = crate::fields::Fp117Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec214 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp214::Fp;
    pub type Fq = crate::fields::Fp214Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec509 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp509::Fp;
    pub type Fq = crate::fields::Fp509Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec1757 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp1757::Fp;
    pub type Fq = crate::fields::Fp1757Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}
