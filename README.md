# hts_utils

C++ library and recipes based on htslib (https://htslib.org) for handling
high-throughput sequencing data.

## Included scripts:

### merge_sams

merge two sam files. If a record for a read exists in both SAM files, report
the record with the "better" score.

```
./merge_sams sam1.[sb]am sam2.[sb]am > out.sam
```

## score_sam

given a truth alignment SAM and a test SAM, check if each record matches the
corresponding record in the truth set.

```
./score_sam truth.[sb]am test.[sb]am

READ	TRUE_REF	TRUE_POS	TEST_REF	TEST_POS	CORRECT?
```
