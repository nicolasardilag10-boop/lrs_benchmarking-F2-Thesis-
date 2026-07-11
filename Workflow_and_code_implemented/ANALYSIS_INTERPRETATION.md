# Analysis Interpretation Guide

## Error-rate interpretation

Lower values indicate fewer alignment errors.

A suitable summary is:

> ONT reads showed consistently higher alignment error rates than PacBio reads across the tested configurations. ONT error rates changed more strongly between configurations, whereas PacBio error rates were lower and more stable.

## Avoid overclaiming significance

Do not say:

> Every individual sample showed a significant difference.

The paired t-test evaluates all paired samples together. It does not calculate a separate significance test for each sample.

## Small sample size

Each technology has only three paired samples:

```text
HG002
HG003
HG004
```

Therefore:

- statistical power is limited,
- p-values should be interpreted cautiously,
- `p = 0.05` is borderline,
- multiple comparisons should be acknowledged.

## Nominal significance

When several comparisons are made without correction, describe the results as **nominal p-values**.

## Mapped-bases interpretation

Absolute mapped bases are influenced by the total number of input bases.

A safe statement is:

> ONT produced more CIGAR-derived mapped bases for HG002, whereas PacBio produced more for HG003 and HG004.

Do not conclude that one technology performed better from absolute mapped-base counts alone.

A stronger comparison combines:

- absolute mapped bases,
- mapped-bases percentage,
- mapped-read percentage,
- alignment error rate.

## Technology-matched comparisons

The most direct comparisons are:

```text
ONT:
mm2-ont versus pbmm2-ont

PacBio:
mm2-pb versus pbmm2-pb
```

## Configuration terminology

A configuration means:

```text
aligner + preset/settings
```

The label `pbmm2-ont` is a project label for the tested pbmm2 configuration applied to ONT input. It should not be presented as an official ONT preset unless documented as such.