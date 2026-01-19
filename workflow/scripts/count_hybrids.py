# Create a pivot table from aradopsis strain sequences and aggregate hybrid ratio using duckdb
from snakemake.script import snakemake
import duckdb

config = snakemake.config

def calc_h_ratio(data_path: str, out_path: str) -> None: # makes new csv file
    """Creates csv file with 'h_ratio' column. 

    First, it will clean up un-necesarry column/rows.
    Second, it will unpivot the chromosome+'_'+position row into a 'chromosome' and 'position' columns.
    Third, it will pivot the genotype values into 'H' (hybrid), 'A' (non-mutant aradposis), 'NA' (not-a-valid-number)

    Args:
        data_path: Input csv file.
        out_path: Output csv file with 'h_ratio' column.

    Returns:
        'None', just saves new csv file.
    """
    
    # create db in memory
    con = duckdb.connect(database=':memory:')

    # Save columns that follow the 'chr_pos' pattern
    initial_rel = con.read_csv(data_path)
    all_cols = initial_rel.columns
    data_cols = [c for c in all_cols if '_' in c and c.split('_')[-1].isdigit()]

    # turn col_list into correct duckdb sql string format
    col_list_sql = ", ".join([f"'{c}'" for c in data_cols])
    
    # remove excluded chromosome positions (centromeres)
    exclude_list = [str(p) for p in config.get("exclude", [])]
    # turn exclude_list into correct duckdb sql string format
    exclude_list_sql = ", ".join([f"'{c}'" for c in exclude_list])
    #data_cols = [c for c in data_cols if c not in set(exclude_list)]
    
    
    # Unpivots position row, then pivots genotype categories, to aggregate H to A ratio
    query = r"""
    -- Remove chr and pos rows (2 and 3)
    WITH cleaned_rel AS (
        SELECT * FROM read_csv('{data_path}')
        WHERE Sample IS NOT NULL AND Sample != ''
    ),
    -- Table with every bp (A/H/NA) having its own row (cols: Sample, chr_pos, genotype)
    long_data AS (
        SELECT 
            Sample,
            header,
            val AS genotype
        FROM cleaned_rel
        UNPIVOT (
            val FOR header IN ({col_list})
        )
    ),
    -- Pivot genotype categories
    pivoted_data as (
        PIVOT long_data 
        ON genotype IN ('H', 'A', 'NA')
        GROUP BY header
    )
    -- Final calculation, split headers here
    SELECT 
        split_part(header, '_', 1)::INTEGER AS chr,
        split_part(header, '_', 2)::INTEGER AS pos,
        H, A, NA,
        (H + A) AS valid_count,
        CASE 
            WHEN header IN ({exclude_list}) THEN NULL
            ELSE ROUND(H::FLOAT / NULLIF((H + A), 0), 3) 
        END AS h_ratio
        --ROUND(H::FLOAT / NULLIF((H + A), 0), 3) AS h_ratio
    FROM pivoted_data
    ORDER BY chr, pos
    """.format(data_path=data_path, col_list=col_list_sql, exclude_list=exclude_list_sql)

    con.sql(query).write_csv(out_path)

if __name__ == "__main__":
    calc_h_ratio(snakemake.input[0], snakemake.output[0])

