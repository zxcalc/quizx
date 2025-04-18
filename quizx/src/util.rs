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
