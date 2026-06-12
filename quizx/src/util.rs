/// The same as max(iterable), but the elements are allowed to be PartialOrd.
///
/// If the iterator contains incomparable items, it will prefer the item that
/// occurs earlier.
pub fn pmax<I>(iterable: I) -> Option<I::Item>
where
    I: IntoIterator,
    I::Item: PartialOrd,
{
    iterable.into_iter().fold(None, {
        |m, it| match m {
            None => Some(it),
            Some(n) => {
                if n < it {
                    Some(it)
                } else {
                    Some(n)
                }
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pmax_integers() {
        assert_eq!(pmax([1, 3, 2]), Some(3));
        assert_eq!(pmax([5, 1, 2, 4]), Some(5));
        assert_eq!(pmax([1]), Some(1));
    }

    #[test]
    fn pmax_empty() {
        assert_eq!(pmax(Vec::<i32>::new()), None);
    }

    #[test]
    fn pmax_floats_with_nan() {
        // NaN is incomparable; pmax prefers earlier element when incomparable.
        let vals = [1.0, f64::NAN, 2.0];
        let result = pmax(vals);
        // NAN compared to 1.0: 1.0 < NAN is false, so keeps 1.0.
        // 1.0 compared to 2.0: 1.0 < 2.0 is true, so takes 2.0.
        assert_eq!(result, Some(2.0));

        // When NaN comes first, it stays (nothing is greater).
        let vals = [f64::NAN, 1.0, 2.0];
        let result = pmax(vals);
        // NAN < 1.0 is false, so keeps NAN.
        // NAN < 2.0 is false, so keeps NAN.
        assert!(result.unwrap().is_nan());
    }

    #[test]
    fn pmax_negative() {
        assert_eq!(pmax([-5, -1, -3]), Some(-1));
    }
}
