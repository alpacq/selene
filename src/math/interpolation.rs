/// Returns the interpolation indexes for a given value and multiplier
///
/// # Arguments
///
/// * `lut_size` - The size of the lookup table
/// * `shift` - The shift value for the lookup table
/// * `value` - The value to interpolate
/// * `multiplier` - The multiplier for the value
///
/// # Returns
///
/// A tuple of the interpolation indexes
pub fn get_interpolation_indexes_and_remainder(
    lut_size: usize,
    shift: i64,
    value: f64,
    multiplier: f64,
) -> (usize, usize, f64) {
    let clamp_bounds = (
        (1.0 - (shift as f64)),
        ((lut_size - 2) as f64) - (shift as f64),
    );
    let index_float = multiplier * value;
    let index_int = index_float.round().clamp(clamp_bounds.0, clamp_bounds.1) as i64;
    let remainder = index_float - index_int as f64;
    let index_int_adjacent = if remainder >= 0.0 {
        index_int + 1
    } else {
        index_int - 1
    };
    (
        (index_int + shift) as usize,
        (index_int_adjacent + shift) as usize,
        remainder,
    )
}

/// Returns the single dimension linear interpolation based on a lookup table
///
/// # Arguments
///
/// * `alpha` - the angle of attack in radians
/// * `lut` - the lookup table of aerodynamic derivatives
///
/// # Returns
///
/// The interpolated value
pub fn lut_interpolation_1d(parameter: f64, lut: &[f64], shift: i64, multiplier: f64) -> f64 {
    let (index, index_adjacent, remainder) =
        get_interpolation_indexes_and_remainder(lut.len(), shift, parameter, multiplier);
    lut[index] + remainder.abs() * (lut[index_adjacent] - lut[index])
}

pub fn lut_interpolation_2d(
    parameter_a: f64,
    parameter_b: f64,
    shift_a: i64,
    shift_b: i64,
    multiplier_a: f64,
    multiplier_b: f64,
    lut: &[&[f64]],
) -> f64 {
    let (a_ind, a_adj, a_rem) =
        get_interpolation_indexes_and_remainder(lut.len(), shift_a, parameter_a, multiplier_a);
    let (b_ind, b_adj, b_rem) =
        get_interpolation_indexes_and_remainder(lut[0].len(), shift_b, parameter_b, multiplier_b);
    let v = lut[a_ind][b_ind] + a_rem.abs() * (lut[a_adj][b_ind] - lut[a_ind][b_ind]);
    let w = lut[a_ind][b_adj] + a_rem.abs() * (lut[a_adj][b_adj] - lut[a_ind][b_adj]);
    v + (w - v) * b_rem.abs()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn positive_remainder_interpolation_indexes() {
        let (index, index_adjacent, _) = get_interpolation_indexes_and_remainder(9, 0, 7.0, 0.2);
        assert_eq!(index, 1);
        assert_eq!(index_adjacent, 2);
    }

    #[test]
    fn clamp_and_positive_remainder_interpolation_indexes() {
        let (index, index_adjacent, _) = get_interpolation_indexes_and_remainder(9, 0, 45.0, 0.2);
        assert_eq!(index, 7);
        assert_eq!(index_adjacent, 8);
    }

    #[test]
    fn negative_remainder_interpolation_indexes() {
        let (index, index_adjacent, _) = get_interpolation_indexes_and_remainder(9, 0, 3.0, 0.2);
        assert_eq!(index, 1);
        assert_eq!(index_adjacent, 0);
    }

    #[test]
    fn clamp_and_negative_remainder_interpolation_indexes() {
        let (index, index_adjacent, _) = get_interpolation_indexes_and_remainder(9, 0, 1.0, 0.2);
        assert_eq!(index, 1);
        assert_eq!(index_adjacent, 0);
    }

    #[test]
    fn negative_value_range() {
        let (index, index_adjacent, _) = get_interpolation_indexes_and_remainder(12, 2, -3.0, 0.2);
        assert_eq!(index, 1);
        assert_eq!(index_adjacent, 2);
    }

    #[test]
    fn negative_value_range_2() {
        let (index, index_adjacent, _) = get_interpolation_indexes_and_remainder(12, 2, -6.0, 0.2);
        assert_eq!(index, 1);
        assert_eq!(index_adjacent, 0);
    }
}
