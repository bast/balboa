use std::cmp::Ordering;
use std::collections::HashMap;

use crate::limits;

pub fn transform_to_spherical(
    num_points: usize,
    aos_c: &[f64],
    l: usize,
    c_to_s_matrix: &[(usize, usize, f64)],
) -> Vec<f64> {
    let spherical_deg = 2 * l + 1;
    let mut aos_s = vec![0.0; spherical_deg * num_points];

    for (i, j, x) in c_to_s_matrix {
        for ipoint in 0..num_points {
            aos_s[i * num_points + ipoint] += x * aos_c[j * num_points + ipoint];
        }
    }

    aos_s
}

pub fn cartesian_to_spherical_matrices() -> HashMap<usize, Vec<(usize, usize, f64)>> {
    let mut map = HashMap::new();

    for l in 0..=limits::MAX_L_VALUE {
        map.insert(l, cartesian_to_spherical_coef(l));
    }

    map
}

// we return f64 to prevent overflow
fn fac(n: usize) -> f64 {
    (0..n).fold(1.0, |acc, i| acc * (i + 1) as f64)
}

// f64 to prevent overflow
fn fac2(n: i32) -> f64 {
    match n.cmp(&0) {
        Ordering::Equal => 1.0,
        Ordering::Greater => {
            let mut r = n as f64;
            let mut i: i32 = (n - 2) as i32;
            while i > 0 {
                r *= i as f64;
                i -= 2;
            }
            r
        }
        Ordering::Less => {
            let mut r = (n + 2) as f64;
            let mut i = n + 4;
            while i < 2 {
                r *= i as f64;
                i += 2;
            }
            1.0 / r
        }
    }
}

// f64 to prevent overflow
fn binom(n: usize, k: usize) -> f64 {
    if k == 0 {
        return 1.0;
    }

    if n == 0 {
        return 0.0;
    }

    let mut m = n as f64;
    let mut b = 1.0;

    for j in 0..k {
        b *= m;
        b /= (j + 1) as f64;
        m -= 1.0;
    }

    b
}

fn cartesian_to_spherical_coef(l: usize) -> Vec<(usize, usize, f64)> {
    if l == 0 {
        return vec![(0, 0, 1.0)];
    }

    if l == 1 {
        return vec![(0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0)];
    }

    let nc = (l + 1) * (l + 2) / 2;
    let ns = 2 * l + 1;

    let mut legendre_coef = vec![0.0; l + 1];
    let mut k = 0;
    while k <= l / 2 {
        let prefactor = (-1.0_f64).powi(k as i32) / (2.0_f64).powi(l as i32);
        let x = binom(l, k) * binom(2 * (l - k), l);
        legendre_coef[l - 2 * k] = prefactor * x;
        k += 1;
    }

    let mut cossin_coef = vec![0.0; (l + 1) * (l + 1)];
    let mut tmat = vec![0.0; nc * ns];
    for m in 0..=l {
        cossin_coef[m] = 1.0;
        for k in 1..=m {
            cossin_coef[k * (l + 1) + m] +=
                cossin_coef[(k - 1) * (l + 1) + m - 1] * (-1.0_f64).powi((k - 1) as i32);
            if m > k {
                cossin_coef[k * (l + 1) + m] += cossin_coef[k * (l + 1) + m - 1];
            }
        }
    }

    for m in 0..=l {
        let mut cm = match m {
            0 => 1.0,
            _ => (2.0 * fac(l - m) / fac(l + m)).sqrt(),
        };
        cm /= (fac2((2 * l - 1) as i32)).sqrt();

        k = (l - m) % 2;
        while k <= (l - m) {
            if m > 0 {
                legendre_coef[k] = legendre_coef[k + 1] * ((k + 1) as f64);
            }
            let cmk = cm * legendre_coef[k];
            let mut i = 0;
            while i <= (l - k - m) / 2 {
                let cmki = cmk * binom((l - k - m) / 2, i);
                for j in 0..=i {
                    let cmkij = cmki * binom(i, j);
                    for n in 0..=m {
                        let mut ix = l - 2 * j - m + n;
                        ix = ix * (ix + 1) / 2 + l + 1 - m - 2 * i;

                        let ilm = match n % 2 {
                            1 => 1 + l - m,
                            _ => 1 + l + m,
                        };

                        tmat[(ilm - 1) * nc + ix - 1] += cmkij * cossin_coef[n * (l + 1) + m];
                    }
                }
                i += 1;
            }
            k += 2;
        }
    }

    let mut tc = Vec::new();
    for i in 0..nc {
        for j in 0..ns {
            let x = tmat[j * nc + i];
            if x.abs() > std::f64::EPSILON {
                tc.push((j, i, x));
            }
        }
    }
    tc
}

#[cfg(test)]
fn are_same(v1: &[(usize, usize, f64)], v2: &[(usize, usize, f64)]) -> bool {
    if v1.len() != v2.len() {
        return false;
    }

    for ((i1, j1, x1), (i2, j2, x2)) in v1.iter().zip(v2.iter()) {
        if i1 != i2 {
            return false;
        }
        if j1 != j2 {
            return false;
        }
        if (x1 - x2).abs() > std::f64::EPSILON {
            return false;
        }
    }

    true
}

#[test]
fn test_l_0() {
    let reference = vec![(0, 0, 1.0)];
    assert!(are_same(&reference, &cartesian_to_spherical_coef(0)));
}

#[test]
fn test_l_1() {
    let reference = vec![(0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0)];
    assert!(are_same(&reference, &cartesian_to_spherical_coef(1)));
}

#[test]
fn test_l_2() {
    let reference = vec![
        (2, 0, -0.2886751345948129),
        (4, 0, 0.5),
        (0, 1, 1.0),
        (3, 2, 1.0),
        (2, 3, -0.2886751345948129),
        (4, 3, -0.5),
        (1, 4, 1.0),
        (2, 5, 0.577350269189626),
    ];
    assert!(are_same(&reference, &cartesian_to_spherical_coef(2)));
}

#[test]
fn test_l_3() {
    let reference = vec![
        (4, 0, -0.15811388300841897),
        (6, 0, 0.20412414523193148),
        (0, 1, 0.6123724356957945),
        (2, 1, -0.15811388300841897),
        (3, 2, -0.38729833462074165),
        (5, 2, 0.5),
        (4, 3, -0.15811388300841897),
        (6, 3, -0.6123724356957945),
        (1, 4, 1.0),
        (4, 5, 0.6324555320336758),
        (0, 6, -0.20412414523193148),
        (2, 6, -0.15811388300841897),
        (3, 7, -0.38729833462074165),
        (5, 7, -0.5),
        (2, 8, 0.6324555320336758),
        (3, 9, 0.25819888974716104),
    ];
    assert!(are_same(&reference, &cartesian_to_spherical_coef(3)));
}

#[test]
fn test_l_4() {
    let reference = vec![
        (4, 0, 0.036596252735569997),
        (6, 0, -0.0545544725589981),
        (8, 0, 0.07216878364870323),
        (0, 1, 0.2886751345948129),
        (2, 1, -0.1091089451179962),
        (5, 2, -0.23145502494313788),
        (7, 2, 0.2041241452319315),
        (4, 3, 0.07319250547113999),
        (8, 3, -0.4330127018922194),
        (1, 4, 0.6123724356957945),
        (3, 4, -0.23145502494313788),
        (4, 5, -0.29277002188456),
        (6, 5, 0.3273268353539886),
        (0, 6, -0.2886751345948129),
        (2, 6, -0.1091089451179962),
        (5, 7, -0.23145502494313788),
        (7, 7, -0.6123724356957945),
        (2, 8, 0.6546536707079772),
        (5, 9, 0.3086066999241839),
        (4, 10, 0.036596252735569997),
        (6, 10, 0.0545544725589981),
        (8, 10, 0.07216878364870323),
        (1, 11, -0.2041241452319315),
        (3, 11, -0.23145502494313788),
        (4, 12, -0.29277002188456),
        (6, 12, -0.3273268353539886),
        (3, 13, 0.3086066999241839),
        (4, 14, 0.09759000729485329),
    ];
    assert!(are_same(&reference, &cartesian_to_spherical_coef(4)));
}
