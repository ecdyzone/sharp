# %%
import pandas as pd
import pyarrow.parquet as pq

input_file = "../data/interim/embeddings.parquet"

# %%
# Load parquet file with pandas
df = pd.read_parquet(input_file)

print(df.head())
print(df.shape)
print(df.columns)
print(df.dtypes)

# %%
# PyArrow low-level inspection (metadata only, fast)
pf = pq.ParquetFile(input_file)

print("num rows:", pf.metadata.num_rows)
print("schema:\n", pf.schema)

# %%
# Load full table with PyArrow (optional sanity check)
table = pq.read_table(input_file)

print(table.slice(0, 5).to_pandas())

# %%
# Basic integrity check
print(f"Length is: {len(df)} lines")
assert len(df) > 0
print("Parquet file looks OK ✔")

