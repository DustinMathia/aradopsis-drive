# Create a pivot table from aradopsis strain sequences and aggregate hybrid ratio using duckdb

import duckdb


relation = duckdb.read_csv(test_file, all_varchar=True)


# Unpivots position row, then pivots genotype categories, to aggregate H to A ratio
query = r"""

-- Remove chr and pos rows (2 and 3)
WITH cleaned_rel AS (
    SELECT * FROM relation
    WHERE Sample IS NOT NULL AND Sample != ''
),

-- Table with every bp (A/H/NA) having its own row (cols: Sample, chr, pos, genotype)
long_data AS (
    SELECT 
        Sample,
        split_part(header, '_', 1)::INTEGER AS chr,
        split_part(header, '_', 2)::INTEGER AS pos,
        val AS genotype
    FROM cleaned_rel
    UNPIVOT (
        val FOR header IN (COLUMNS('^.+_+\d'))
    )
),

-- Split genotype column into three category columns: H, A, NA
pivoted_data as (
    PIVOT long_data 
    ON genotype IN ('H', 'A', 'NA')
    GROUP BY chr, pos
)

-- make two new columns, total bp's that aren't NA at that pos, and ratio of H vs. A at that pos
SELECT 
    *,
    (H + A) AS valid_count,
    ROUND(H::FLOAT / NULLIF((H + A), 0), 3) AS h_ratio
FROM pivoted_data
ORDER BY chr, pos
"""

# Execute query
tidy_data = duckdb.sql(query)