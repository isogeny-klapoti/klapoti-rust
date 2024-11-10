use crate::{linalg::matrix::Matrix, util::Big};
use rug::Integer;

use super::lattice::Lattice;

#[derive(Clone, Debug, PartialEq)]
pub struct QuaternionOrder {
    pub lattice: Lattice,
}

impl QuaternionOrder {
    pub fn new(lattice: Lattice) -> Self {
        QuaternionOrder { lattice }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ExtremalMaximalOrder {
    pub order: QuaternionOrder,
    pub q: Integer, // absolute value of sqrt of i
}

impl ExtremalMaximalOrder {
    pub fn new(order: QuaternionOrder, q: Integer) -> Self {
        ExtremalMaximalOrder { order, q }
    }
}

pub fn standard_extremal_order() -> ExtremalMaximalOrder {
    let basis = Matrix::zeros(4, 4);
    let mut lat = Lattice::new(basis.clone(), 2.big());
    lat.basis[(0, 0)] = 2.big();
    lat.basis[(0, 3)] = 1.big();
    lat.basis[(1, 1)] = 2.big();
    lat.basis[(1, 2)] = 1.big();
    lat.basis[(2, 2)] = 1.big();
    lat.basis[(3, 3)] = 1.big();
    let order = QuaternionOrder::new(lat);

    ExtremalMaximalOrder::new(order, 1.big())
}

pub fn standard_extremal_from(order: QuaternionOrder) -> ExtremalMaximalOrder {
    ExtremalMaximalOrder::new(order, 1.big())
}
