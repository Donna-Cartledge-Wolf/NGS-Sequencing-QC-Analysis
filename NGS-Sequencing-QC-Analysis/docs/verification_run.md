# Verification record — v1.0.0

The project was executed against the included synthetic FASTQ dataset and
all 13 standard-library unit tests passed in the build environment.

Expected reference checkpoints (verified by automated tests):

| Checkpoint | Value |
|---|---:|
| Reads | 12 |
| Total bases | 720 |
| Q30+ bases | 611 |
| Reads with educational review flags | 6 |
| Bases removed by terminal Q20 trimming | 14 |
| SIM_02 expected errors before trimming | 1.21517893 |
| SIM_02 expected errors after trimming | 0.01517893 |
| SIM_03 expected errors | 0.06000000 |

The 14 removed bases include SIM_02's 12-base Q10 tail and one low-quality
terminal base from each of SIM_11 and SIM_12.

Run `python -m unittest discover -s tests -v` locally to reproduce results.
These tests verify software behavior on one synthetic dataset and edge cases.
They are not assay validation, certification, or performance qualification.
