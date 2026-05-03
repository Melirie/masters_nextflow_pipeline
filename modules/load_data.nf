process load_data {
    input:
    val url
    val filename

    output:
    path "${filename}"

    script:
    """
    python3 -c "
import pooch
# Everything below must be at the very start of the line (no spaces!)
pooch.retrieve(
    url='${url}',
    known_hash=None,
    fname='${filename}',
    path='.',
    progressbar=True
)
"
    """
}
