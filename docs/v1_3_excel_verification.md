# Version 1.3 Excel verification exercise

Use the Python-generated `results/synthetic_gapped_alignments.csv` as the
reference output. Save a separate `.xlsx` working copy if you want to add
formulas and worksheets. Do not edit the committed Python output CSV directly.

## Exercise A: insertion in the query (global)

Open the row named `INSERTION_IN_QUERY` with `method=global`.

Expected alignment from the **synthetic** example:

```text
Reference: ACGT-TGCA
Query:     ACGTATGCA
Operation: 4=1I4=
```

There are eight matching alignment columns and one insertion column.

| Metric | Excel formula | Expected |
|---|---|---:|
| Aligned columns | `=LEN("ACGT-TGCA")` | 9 |
| Matching columns | `=8` | 8 |
| Percent identity (including gaps) | `=8/9*100` | 88.88888889 |
| Default alignment score | `=8*2+1*(-2)` | 14 |

## Exercise B: deletion in the query (global)

```text
Reference: ACGTTGCA
Query:     ACG-TGCA
Operation: 3=1D4=
```

| Metric | Excel formula | Expected |
|---|---|---:|
| Aligned columns | `=LEN("ACGTTGCA")` | 8 |
| Matching columns | `=7` | 7 |
| Percent identity (including gaps) | `=7/8*100` | 87.5 |
| Default alignment score | `=7*2+1*(-2)` | 12 |

**Reminder:** a gap in the query is a deletion *relative to the reference*.
The alignment does not establish that a biological deletion occurred.

## Exercise C: why local and global alignment differ

Compare the `LOCAL_FLANKS` rows in `synthetic_gapped_alignments.csv`:

- Global alignment score: 6. Its selected alignment covers both complete
  strings and includes mismatches and gaps.
- Local alignment score: 16. It retains an exact eight-base shared region,
  `ACCTGATC`, and omits nonmatching flanks.

A 100% identity for the local span is **not** 100% identity for the whole query.
Also inspect the reference and query coordinates (1-based inclusive).

## Quality review questions

1. Does a nine-column alignment with eight matches prove a true insertion?
2. Why might homopolymer gaps have multiple equally optimal placements?
3. Would a changed scoring scheme necessarily change biological truth?
4. Why does the identity denominator need to be specified in a report?

All calculations and interpretations are educational examples. They are
not independent biological confirmation or a validated sequencing assay.
