#[derive(Clone)]
pub struct Shell {
    pub coordinates: (f64, f64, f64),
    pub l: usize,
    pub primitive_exponents: Vec<f64>,
    pub contraction_coefficients: Vec<f64>,
}

pub struct Basis {
    pub num_ao_cartesian: usize,
    pub num_ao_spherical: usize,
    pub shells: Vec<Shell>,
    pub shell_extents_squared: Vec<f64>,
}

impl Basis {
    pub fn new(shells: &[Shell]) -> Basis {
        let num_ao_cartesian = shells
            .iter()
            .fold(0, |sum, s| sum + ((s.l + 1) * (s.l + 2) / 2));
        let num_ao_spherical = shells.iter().fold(0, |sum, s| sum + (2 * s.l + 1));

        Basis {
            num_ao_cartesian,
            num_ao_spherical,
            shells: shells.to_vec(),
            shell_extents_squared: shells.iter().map(shell_extent_squared).collect(),
        }
    }
}

// threshold and factors match Dalton implementation, see also pink book
fn shell_extent_squared(shell: &Shell) -> f64 {
    let f = vec![1.0, 1.3333, 1.6, 1.83, 2.03, 2.22, 2.39, 2.55, 2.70, 2.84];
    let shell_screening_threshold: f64 = 2.0e-12;

    let mut r: f64 = 0.0;
    for (e, c) in shell
        .primitive_exponents
        .iter()
        .zip(shell.contraction_coefficients.iter())
    {
        r = r.max((c.abs().ln() - shell_screening_threshold.ln()) / e);
    }

    if shell.l < 10 {
        r = r.sqrt() * f[shell.l];
    } else {
        r = 1.0e10;
    }

    r * r
}
