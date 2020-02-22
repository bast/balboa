pub fn get_ijk_list(m: usize) -> Vec<(usize, usize, usize)> {
    let mut l = Vec::new();

    for a in 1..(m + 2) {
        for b in 1..(a + 1) {
            l.push((m + 1 - a, a - b, b - 1));
        }
    }

    return l;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ijk() {
        let reference = vec![
            (3, 0, 0),
            (2, 1, 0),
            (2, 0, 1),
            (1, 2, 0),
            (1, 1, 1),
            (1, 0, 2),
            (0, 3, 0),
            (0, 2, 1),
            (0, 1, 2),
            (0, 0, 3),
        ];
        assert_eq!(reference, get_ijk_list(3));
    }
}
