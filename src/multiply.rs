#![allow(clippy::many_single_char_names)]

use crate::diff;
use crate::limits;

pub fn multiply_batch(
    cartesian_orders: (usize, usize, usize),
    geo_derv_orders: (usize, usize, usize),
    gaussians: &[Vec<f64>],
    pxs: &[f64],
    pys: &[f64],
    pzs: &[f64],
    aos_c: &mut [f64],
) {
    let map = diff::differentiate(cartesian_orders, geo_derv_orders);

    for (key, value) in map {
        let i = key[0] as i32;
        let j = key[1] as i32;
        let k = key[2] as i32;
        let p = key[3];

        let v = value as f64;

        for ipoint in 0..limits::BATCH_LENGTH {
            aos_c[ipoint] += v
                * gaussians[ipoint][p]
                * pxs[ipoint].powi(i)
                * pys[ipoint].powi(j)
                * pzs[ipoint].powi(k);
        }
    }
}
