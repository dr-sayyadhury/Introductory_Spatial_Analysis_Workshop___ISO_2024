import pandas as pd

# Description: Read parquet file and return a pandas dataframe
def check_parquet(filename):
    df = pd.read_parquet(filename)
    for col in df.columns:
        if df[col].dtype == 'object':  # Check if column is of object type
            df[col] = df[col].apply(lambda x: x.decode('utf-8') if isinstance(x, bytes) else x)
    return df

