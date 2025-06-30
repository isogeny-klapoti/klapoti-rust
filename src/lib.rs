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

pub mod ec64 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp64::Fp;
    pub type Fq = crate::fields::Fp64Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec128 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp128::Fp;
    pub type Fq = crate::fields::Fp128Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec256 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp256::Fp;
    pub type Fq = crate::fields::Fp256Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec512 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp512::Fp;
    pub type Fq = crate::fields::Fp512Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec768 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp768::Fp;
    pub type Fq = crate::fields::Fp768Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec1024 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp1024::Fp;
    pub type Fq = crate::fields::Fp1024Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec1536 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp1536::Fp;
    pub type Fq = crate::fields::Fp1536Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}

pub mod ec2048 {
    use crate::util::Big;
    pub type Fp = crate::fields::Fp2048::Fp;
    pub type Fq = crate::fields::Fp2048Ext::Fp2;
    crate::ec::eccore::define_ec_core! {}
    crate::theta::theta::define_theta_structure! {}
    crate::schemes::klapoti::define_klapoti! {}
}