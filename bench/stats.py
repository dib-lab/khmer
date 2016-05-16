import pandas as pd
import sarparse

def runtime(df):
    return df.iloc[-1]['time'] - df.iloc[0]['time']

df_nosh = pd.DataFrame(sarparse.parse_pidstat("output_nosh.pidstat.gz"))
df_nosh['time'] = pd.to_datetime(df_nosh['time'])
df_nosh.set_index("time")

df_sh = pd.DataFrame(sarparse.parse_pidstat("output_sh.pidstat.gz"))
df_sh['time'] = pd.to_datetime(df_sh['time'])
df_sh.set_index("time")
