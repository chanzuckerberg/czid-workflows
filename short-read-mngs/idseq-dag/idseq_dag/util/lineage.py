from typing import List, Sequence, Iterator

INVALID_CALL_BASE_ID = -(10**8)
NULL_SPECIES_ID = "-100"
NULL_GENUS_ID = "-200"
NULL_FAMILY_ID = "-300"
NULL_LINEAGE = (NULL_SPECIES_ID, NULL_GENUS_ID, NULL_FAMILY_ID)

DEFAULT_BLACKLIST_S3 = 's3://czid-public-references/taxonomy/2018-04-01-utc-1522569777-unixtime__2018-04-04-utc-1522862260-unixtime/taxon_blacklist.txt'
DEFAULT_WHITELIST_S3 = 's3://czid-public-references/taxonomy/2020-02-10/respiratory_taxon_whitelist.txt'

# Notes for lineage and taxonomy ID functions:
#
# A hit will normally have (positive) NCBI taxon IDs at all levels of the
# hierarchy, but we use an artificial negative taxon ID if we have determined
#  that the alignment is not specific at the taxonomy level under
# consideration. This happens when a read's multiple reference matches do not
#  agree on taxon ID at the given level.
#
# For example, a read may match 5 references that all belong to different
# species (e.g. Escherichia albertii, Escherichia vulneris, Escherichia coli,
#  ...), but to the same genus (Escherichia). In this case, we use the taxon
# ID for the genus (Escherichia) at the genus-level, but we populate the
# species-level with an artificial negative ID. The artificial ID is defined
# based on a negative base ( INVALID_CALL_BASE_ID), the taxon level (e.g. 2
# for genus), and the valid parent ID (e.g. genus Escherichia's taxon ID).
#
# See also comments under call_hits_m8().


def _cleaned_taxid_lineage(taxid_lineage: Sequence[str], hit_taxid: str, hit_level: int) -> Iterator[str]:
    """Take the taxon lineage and mark meaningless calls with fake taxids."""
    # This assumption is being made in postprocessing
    assert len(taxid_lineage) == 3
    for tax_level, taxid in enumerate(taxid_lineage, 1):
        if tax_level >= hit_level:
            yield taxid
        else:
            yield tax_level * INVALID_CALL_BASE_ID - int(hit_taxid)


def _fill_missing_calls(cleaned_lineage: Sequence[str]) -> List[str]:
    """Replace missing calls with virtual taxids as shown in
    fill_missing_calls_tests. Replaces the negative call IDs with a
    calculated negative value that embeds the next higher positive call in it
    to indicate (1) the call is non-specific (negative/artificial) and (2)
    there is a positive specific call above it on the taxonomy tree.

    Ex: with species_id, genus_id, family_id
    (55, -200, 1534) => (55, -200001534, 1534)
    (-100, 5888, -300) => (-100005888, 5888, -300)
    """
    result = list(cleaned_lineage)
    tax_level = len(cleaned_lineage)
    closest_real_hit_just_above_me = -1

    for tax_level in range(len(cleaned_lineage), 0, -1):
        me = int(cleaned_lineage[tax_level - 1])
        if me >= 0:
            closest_real_hit_just_above_me = me
        elif closest_real_hit_just_above_me >= 0 and _blank(me):
            result[tax_level - 1] = str(tax_level * INVALID_CALL_BASE_ID - closest_real_hit_just_above_me)
    return result


def _blank(taxid_int: int):
    return 0 > taxid_int > INVALID_CALL_BASE_ID


def fill_missing_calls_tests():
    # TODO: Move these to a test file later
    # -200 => -200001534
    assert _fill_missing_calls(("55", "-200", "1534")) == ["55", "-200001534", "1534"]
    # -100 => -100005888
    assert _fill_missing_calls(("-100", "5888", "-300")) == ["-100005888", "5888", "-300"]
    # -100 => -100001534, -200 => -200001534
    assert _fill_missing_calls(("-100", "-200",
                               "1534")) == ["-100001534", "-200001534", "1534"]
    # no change
    assert _fill_missing_calls(("55", "-200", "-300")) == ["55", "-200", "-300"]


def validate_taxid_lineage(taxid_lineage: Sequence[str], hit_taxid: str, hit_level: int) -> List[str]:
    cleaned = _cleaned_taxid_lineage(taxid_lineage, hit_taxid, hit_level)
    return _fill_missing_calls(list(cleaned))
