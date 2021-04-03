use std::collections::HashMap;

pub type DiffMap = HashMap<[usize; 4], usize>;

pub fn differentiate(
    cartesian_orders: (usize, usize, usize),
    geo_diff_orders: (usize, usize, usize),
) -> DiffMap {
    let mut map = HashMap::new();
    map.insert(
        [
            cartesian_orders.0,
            cartesian_orders.1,
            cartesian_orders.2,
            0,
        ],
        1,
    );

    for _ in 1..=geo_diff_orders.0 {
        map = differentiate_direction(map, 0);
    }
    for _ in 1..=geo_diff_orders.1 {
        map = differentiate_direction(map, 1);
    }
    for _ in 1..=geo_diff_orders.2 {
        map = differentiate_direction(map, 2);
    }

    map
}

fn differentiate_direction(ts: DiffMap, direction: usize) -> DiffMap {
    let mut ts_new = HashMap::new();

    for (key, value) in ts {
        let power = key[direction];

        // product rule to every element we find
        // along the direction of differentiation
        if power > 0 {
            let mut t = key;
            t[direction] -= 1;
            let old_value = match ts_new.get(&t) {
                Some(&v) => v,
                _ => 0,
            };
            ts_new.insert(t, old_value + value * power);
        }

        // this terms comes from differentiating the exp part
        let mut t = key;
        t[direction] += 1;
        t[3] += 1;
        let old_value = match ts_new.get(&t) {
            Some(&v) => v,
            _ => 0,
        };
        ts_new.insert(t, old_value + value);
    }

    ts_new
}

#[test]
fn test_s_xxyy() {
    let map = differentiate((0, 0, 0), (2, 2, 0));

    let mut map_reference: DiffMap = HashMap::new();
    map_reference.insert([0, 0, 0, 2], 1);
    map_reference.insert([0, 2, 0, 3], 1);
    map_reference.insert([2, 0, 0, 3], 1);
    map_reference.insert([2, 2, 0, 4], 1);

    assert_eq!(map, map_reference);
}

#[test]
fn test_s_xxxxxx() {
    let map = differentiate((0, 0, 0), (6, 0, 0));

    let mut map_reference: DiffMap = HashMap::new();
    map_reference.insert([4, 0, 0, 5], 15);
    map_reference.insert([2, 0, 0, 4], 45);
    map_reference.insert([6, 0, 0, 6], 1);
    map_reference.insert([0, 0, 0, 3], 15);

    assert_eq!(map, map_reference);
}
