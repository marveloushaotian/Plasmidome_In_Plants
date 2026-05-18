# Master table provenance investigation

This note records the evidence found while tracing how the current canonical
master table was produced.

## Files inspected

Current canonical final table:

`Collect/NCBI_4395_Batch/Master_Table/final/07_contig_annotation_master_table.csv`

Current reproducible intermediate table:

`Collect/NCBI_4395_Batch/Master_Table/intermediate/04_contig_sample_annotations.provirus_overlap.csv`

Older downstream result table still present in `Result`:

`Result/NCBI_4395_Batch/05_Tree/Genus_Level/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv`

MOB-typer contig mobility input:

`Collect/NCBI_4395_Batch/Master_Table/input/05_mobtyper_contig_mobility.raw.csv`

Contig category reference table:

`Collect/NCBI_4395_Batch/Master_Table/input/06_locus_contig_type_reference.raw.tsv`

Network table derived after the current final table:

`Result/NCBI_4395_Batch/07_Network/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected_with_Transferred.csv`

## Observed table chain

The evidence supports this chain:

```text
input/01_contig_sample_map.raw.csv
input/02_defense_annotations.raw.csv
input/03_amr_annotations.raw.csv
input/04_sample_metadata.raw.csv
        |
        v
intermediate/01_defense_annotations.by_contig.csv
intermediate/02_amr_annotations.by_contig.csv
intermediate/03_contig_sample_annotations.merged.csv
intermediate/04_contig_sample_annotations.provirus_overlap.csv
        |
        v
parse GeNomad_Contig_Type from Contig_ID
rename GTDB taxonomy to CRBC taxonomy
join input/05_mobtyper_contig_mobility.raw.csv
join input/06_locus_contig_type_reference.raw.tsv as Locus_Mapped_Contig_Type
drop redundant Contig_Type2
        |
        v
final/07_contig_annotation_master_table.csv
        |
        v
Result/.../Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected_with_Transferred.csv
```

## Step A: current script can reproduce 01-05

`101_build_master_table.py` rebuilds the four normalized intermediate outputs
and the final master table:

- `01_defense_annotations.by_contig.csv`
- `02_amr_annotations.by_contig.csv`
- `03_contig_sample_annotations.merged.csv`
- `04_contig_sample_annotations.provirus_overlap.csv`
- `07_contig_annotation_master_table.csv`

The final table is now generated from the normalized inputs, not copied from an
existing curated final file.

## Step B: 04 to the older 60-column result table

The old 60-column table has the same rows and same row order as 04:

| Table | Rows | Columns |
| --- | ---: | ---: |
| `04_contig_sample_annotations.provirus_overlap.csv` | 176235 | 63 |
| `Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv` | 176235 | 60 |

`Sample_ID + Contig_ID` is unique in both tables, and all 176235 keys match.

The old 60-column table appears to be derived from 04 by:

1. Parsing `Contig_ID` to create `GeNomad_Contig_Type`.
2. Parsing `Contig_ID` to create `Contig_Type2`.
3. Renaming GTDB taxonomy columns to `*_CRBC`.
4. Creating `Genus_CRBC_Updated`.
5. Joining `predicted_mobility` from `05_mobtyper_contig_mobility.raw.csv`.
6. Dropping the original raw taxonomy columns and the original `* (gtdb)` columns.

### Contig type evidence

`Contig_Type2` can be reproduced exactly from `Contig_ID`:

| `Contig_Type2` | Rows |
| --- | ---: |
| Chromosome | 152827 |
| Plasmid | 22424 |
| Virus | 984 |

`GeNomad_Contig_Type` is almost the same, but keeps provirus information:

| `GeNomad_Contig_Type` | Rows |
| --- | ---: |
| Chromosome | 148421 |
| Plasmid | 22376 |
| Chromosome\|Provirus | 4406 |
| Virus | 984 |
| Plasmid\|provirus | 48 |

### GTDB/CRBC taxonomy evidence

The `*_CRBC` columns in the old 60-column table match the `* (gtdb)` columns
from 04. For example:

- `Kingdom_CRBC` matches `Kingdom (gtdb)` for all 176235 rows.
- `Phylum_CRBC` matches `Phylum (gtdb)` for all 176235 rows.
- `Class_CRBC` matches `Class (gtdb)` for all 176235 rows.
- `Order_CRBC` matches `Order (gtdb)` for all 176235 rows.
- `Family_CRBC` matches `Family (gtdb)` for all 176235 rows.
- `Species_CRBC` matches `Species (gtdb)` for all 176235 rows.

`Genus_CRBC` matches `Genus (gtdb)` for all 176235 rows.

### Genus suffix cleanup evidence

`Genus_CRBC_Updated` removes GTDB genus suffixes such as `_A`, `_B`, `_AB`, etc.

Examples:

| `Genus_CRBC` | `Genus_CRBC_Updated` |
| --- | --- |
| `Steroidobacter_A` | `Steroidobacter` |
| `Pseudoxanthomonas_A` | `Pseudoxanthomonas` |
| `Pseudomonas_E` | `Pseudomonas` |
| `Paenibacillus_AB` | `Paenibacillus` |

Empty genus values are converted to `unknown`.

### Mobility evidence

`predicted_mobility` comes from:

`Collect/NCBI_4395_Batch/Master_Table/input/05_mobtyper_contig_mobility.raw.csv`

That file has 152581 rows, all with unique `sample_id`. When it is joined to the
old 60-column table by:

```text
old_table.Contig_ID == 05_mobtyper_contig_mobility.raw.csv sample_id
```

all 152581 non-empty `predicted_mobility` values are reproduced exactly.

This is also supported by the shell history from 2025-11-26:

```text
python -c "... m.merge(t[['sample_id','predicted_mobility']], left_on='Contig_ID', right_on='sample_id', how='left') ..."
```

## Step C: old 60-column table to current 05 table

The current final master table can be reproduced from the normalized 04 table,
`05_mobtyper_contig_mobility.raw.csv`, and
`06_locus_contig_type_reference.raw.tsv`.

| Table | Rows | Columns |
| --- | ---: | ---: |
| Old 60-column result table | 176235 | 60 |
| Current `07_contig_annotation_master_table.csv` | 155363 | 60 |

All 155363 `Sample_ID + Contig_ID` keys in 05 are present in the old 60-column
table. All 155363 `Contig_ID` values in 05 also match `Downstream_ID` values in
`06_locus_contig_type_reference.raw.tsv`.

The exact final transformation is:

1. Read `04_contig_sample_annotations.provirus_overlap.csv`.
2. Parse `GeNomad_Contig_Type` from `Contig_ID`.
3. Rename the GTDB taxonomy columns to the CRBC taxonomy columns.
4. Create `Genus_CRBC_Updated`.
5. Join `predicted_mobility` from `05_mobtyper_contig_mobility.raw.csv`.
6. Read `06_locus_contig_type_reference.raw.tsv`.
7. Keep only rows where:

   ```text
   table.Contig_ID == locus_mapping.Downstream_ID
   ```

8. Add:

   ```text
   Locus_Mapped_Contig_Type = title_case(locus_mapping.category)
   ```

   This converts `chromosome`, `plasmid`, and `ambiguous` to `Chromosome`,
   `Plasmid`, and `Ambiguous`.

9. Drop `Contig_Type2`, because it is only the simplified form of
   `GeNomad_Contig_Type`.

The generated final table matches the previous 05 table after ignoring the
single removed `Contig_Type2` column.

The current no-`Contig_Type2` 05 table has MD5:

```text
b0ed44653691fecad8fe0ababc3f64d9
```

### Removed rows

Rows removed from the old 60-column table before the current 05:

| Removed category | Rows |
| --- | ---: |
| `Contig_Type2 = Chromosome` | 20740 |
| `Contig_Type2 = Plasmid` | 92 |
| `Contig_Type2 = Virus` | 40 |

All provirus-labelled rows were removed:

| `GeNomad_Contig_Type` | Removed rows |
| --- | ---: |
| Chromosome\|Provirus | 4406 |
| Plasmid\|provirus | 48 |

The row-filtering rule is explained by the exact join to
`06_locus_contig_type_reference.raw.tsv`: rows absent from that reference table
are not kept in 05.

### `Locus_Mapped_Contig_Type` evidence

`Locus_Mapped_Contig_Type` in the current 05 table has:

| `Locus_Mapped_Contig_Type` | Rows |
| --- | ---: |
| Chromosome | 132305 |
| Plasmid | 20819 |
| Ambiguous | 2239 |

Relationship between `Contig_Type2` and `Locus_Mapped_Contig_Type`:

| `Contig_Type2` | `Locus_Mapped_Contig_Type` | Rows |
| --- | --- | ---: |
| Chromosome | Chromosome | 131370 |
| Plasmid | Plasmid | 20810 |
| Plasmid | Ambiguous | 1522 |
| Virus | Chromosome | 935 |
| Chromosome | Ambiguous | 717 |
| Virus | Plasmid | 9 |

The rule that assigns `Locus_Mapped_Contig_Type` is the `category` column from
`06_locus_contig_type_reference.raw.tsv`.

## Step D: current 05 to transferred-plasmid network table

The network table:

`Result/NCBI_4395_Batch/07_Network/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected_with_Transferred.csv`

has 155363 rows and 62 columns. It matches the current 05 table plus:

`Transferred_Plasmid`

This table is downstream of 05, not the source of 05.

## Current conclusion

The current reproducible path is:

```text
normalized input files
-> intermediate 01-04
-> parsed GeNomad_Contig_Type
-> CRBC taxonomy columns
-> MOB-typer predicted_mobility
-> input/06_locus_contig_type_reference.raw.tsv category as Locus_Mapped_Contig_Type
-> final master table without Contig_Type2
```

`101_build_master_table.py` is now the deterministic script from the four
normalized input files to the current 05 table.
