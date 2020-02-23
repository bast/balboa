use crate::limits;

// we return f64 to prevent overflow
fn fac(n: usize) -> f64 {
    return (0..n).fold(1.0, |acc, i| acc * (i + 1) as f64);
}

// f64 to prevent overflow
fn fac2(n: i32) -> f64 {
    if n < 0 {
        let mut r = (n + 2) as f64;
        let mut i = n + 4;
        while i < 2 {
            r *= i as f64;
            i += 2;
        }
        return 1.0 / r;
    } else if n == 0 {
        return 1.0;
    } else {
        let mut r = n as f64;
        let mut i: i32 = (n - 2) as i32;
        while i > 0 {
            r *= i as f64;
            i -= 2;
        }
        return r;
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
        b = b * m;
        b = b / ((j + 1) as f64);
        m -= 1.0;
    }

    return b;
}

fn cartesian_to_spherical_coef(l: usize) -> Vec<Vec<f64>> {
    if l == 0 {
        return vec![vec![1.0]];
    }

    if l == 1 {
        return vec![
            vec![0.0, 0.0, 1.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
        ];
    }

    let nc = (l + 1) * (l + 2) / 2;
    let ns = 2 * l + 1;

    let mut legendre_coef = vec![0.0; l + 1];
    let mut k = 0;
    while k <= l / 2 {
        let prefactor = (-1.0 as f64).powi(k as i32) / (2.0 as f64).powi(l as i32);
        let x = binom(l, k) * binom(2 * (l - k), l);
        legendre_coef[l - 2 * k] = prefactor * x;
        k += 1;
    }

    let mut cossin_coef = vec![0.0; (l + 1) * (l + 1)];
    let mut tmat = vec![0.0; nc * ns];
    for m in 0..(l + 1) {
        cossin_coef[m] = 1.0;
        for k in 1..(m + 1) {
            cossin_coef[k * (l + 1) + m] +=
                cossin_coef[(k - 1) * (l + 1) + m - 1] * (-1.0 as f64).powi((k - 1) as i32);
            if m > k {
                cossin_coef[k * (l + 1) + m] += cossin_coef[k * (l + 1) + m - 1];
            }
        }
    }

    for m in 0..(l + 1) {
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
                for j in 0..(i + 1) {
                    let cmkij = cmki * binom(i, j);
                    for n in 0..(m + 1) {
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
        let mut ts = Vec::new();
        for j in 0..ns {
            ts.push(tmat[j * nc + i]);
        }
        tc.push(ts);
    }

    return tc;
}

pub fn cartesian_to_spherical_matrices() -> Vec<Vec<Vec<f64>>> {
    let mut m = Vec::new();

    for l in 0..(limits::MAX_L_VALUE + 1) {
        m.push(cartesian_to_spherical_coef(l));
    }

    return m;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn are_same(vec_vec1: Vec<Vec<f64>>, vec_vec2: Vec<Vec<f64>>) -> bool {
        if vec_vec1.len() != vec_vec2.len() {
            return false;
        }
        for (vec1, vec2) in vec_vec1.iter().zip(vec_vec2.iter()) {
            if vec1.len() != vec2.len() {
                return false;
            }
            for (&e1, &e2) in vec1.iter().zip(vec2.iter()) {
                let d = e1 - e2;
                if d.abs() > 1.0e-8 {
                    return false;
                }
            }
        }
        return true;
    }

    #[test]
    fn test_l_0() {
        let reference = vec![vec![1.0]];
        assert!(are_same(reference, cartesian_to_spherical_coef(0)));
    }

    #[test]
    fn test_l_1() {
        let reference = vec![
            vec![0.0, 0.0, 1.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
        ];
        assert!(are_same(reference, cartesian_to_spherical_coef(1)));
    }

    #[test]
    fn test_l_2() {
        let reference = vec![
            vec![0.0, 0.0, -0.2886751345948129, 0.0, 0.5],
            vec![1.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, -0.2886751345948129, 0.0, -0.5],
            vec![0.0, 1.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.577350269189626, 0.0, 0.0],
        ];
        assert!(are_same(reference, cartesian_to_spherical_coef(2)));
    }

    #[test]
    fn test_l_3() {
        let reference = vec![
            vec![
                0.0,
                0.0,
                0.0,
                0.0,
                -0.15811388300841897,
                0.0,
                0.20412414523193148,
            ],
            vec![
                0.6123724356957945,
                0.0,
                -0.15811388300841897,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            vec![
                0.0,
                0.0,
                0.0,
                -0.38729833462074165,
                0.0,
                0.4999999999999999,
                0.0,
            ],
            vec![
                0.0,
                0.0,
                0.0,
                0.0,
                -0.15811388300841897,
                0.0,
                -0.6123724356957945,
            ],
            vec![0.0, 0.9999999999999998, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.6324555320336758, 0.0, 0.0],
            vec![
                -0.20412414523193148,
                0.0,
                -0.15811388300841897,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            vec![
                0.0,
                0.0,
                0.0,
                -0.38729833462074165,
                0.0,
                -0.4999999999999999,
                0.0,
            ],
            vec![0.0, 0.0, 0.6324555320336758, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.25819888974716104, 0.0, 0.0, 0.0],
        ];
        assert!(are_same(reference, cartesian_to_spherical_coef(3)));
    }

    #[test]
    fn test_l_4() {
        let reference = vec![
            vec![
                0.0,
                0.0,
                0.0,
                0.0,
                0.036596252735569997,
                0.0,
                -0.0545544725589981,
                0.0,
                0.07216878364870323,
            ],
            vec![
                0.2886751345948129,
                0.0,
                -0.1091089451179962,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            vec![
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                -0.23145502494313788,
                0.0,
                0.2041241452319315,
                0.0,
            ],
            vec![
                0.0,
                0.0,
                0.0,
                0.0,
                0.07319250547113999,
                0.0,
                0.0,
                0.0,
                -0.4330127018922194,
            ],
            vec![
                0.0,
                0.6123724356957945,
                0.0,
                -0.23145502494313788,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            vec![
                0.0,
                0.0,
                0.0,
                0.0,
                -0.29277002188456,
                0.0,
                0.3273268353539886,
                0.0,
                0.0,
            ],
            vec![
                -0.2886751345948129,
                0.0,
                -0.1091089451179962,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            vec![
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                -0.23145502494313788,
                0.0,
                -0.6123724356957945,
                0.0,
            ],
            vec![0.0, 0.0, 0.6546536707079772, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.3086066999241839, 0.0, 0.0, 0.0],
            vec![
                0.0,
                0.0,
                0.0,
                0.0,
                0.036596252735569997,
                0.0,
                0.0545544725589981,
                0.0,
                0.07216878364870323,
            ],
            vec![
                0.0,
                -0.2041241452319315,
                0.0,
                -0.23145502494313788,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            vec![
                0.0,
                0.0,
                0.0,
                0.0,
                -0.29277002188456,
                0.0,
                -0.3273268353539886,
                0.0,
                0.0,
            ],
            vec![0.0, 0.0, 0.0, 0.3086066999241839, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.09759000729485329, 0.0, 0.0, 0.0, 0.0],
        ];
        assert!(are_same(reference, cartesian_to_spherical_coef(4)));
    }
}
