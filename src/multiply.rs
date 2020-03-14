use crate::diff;
use crate::limits;

pub fn multiply_batch(
    cartesian_orders: (usize, usize, usize),
    geo_derv_orders: (usize, usize, usize),
    aos_c: &mut Vec<f64>,
    gaussians: &[Vec<f64>],
    pxs: &[f64],
    pys: &[f64],
    pzs: &[f64],
    ao_offset: usize,
    point_offset: usize,
) {
    let map = diff::differentiate(cartesian_orders, geo_derv_orders);

    for (key, value) in map {
        let i = key[0];
        let j = key[1];
        let k = key[2];
        let p = key[3];

        for ipoint in 0..limits::BATCH_LENGTH {
            aos_c[ao_offset + ipoint] += (value as f64)
                * gaussians[point_offset + ipoint][p]
                * pxs[point_offset + ipoint].powi(i as i32)
                * pys[point_offset + ipoint].powi(j as i32)
                * pzs[point_offset + ipoint].powi(k as i32);
        }
    }
}
