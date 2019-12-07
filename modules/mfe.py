import sys
sys.path.append("/usr/local/lib/python3.7/site-packages")
import RNA
import pandas as pd
import dask.dataframe as dd
import numpy as np

def mfe_sliding_window(seq, window_size=100, step_size=50):
  seq_len = len(seq)
  if seq_len < window_size:
    return (np.nan,)*4
  else:
    # Get window start positions
    window_start_positions = list(range(0, seq_len - window_size + 1, step_size))
    if (seq_len - window_size) % step_size: 
      # Slide the window to the 5' end for the last piece
      window_start_positions += [seq_len - window_size] 
    # Calculate and accumulate mfe
    mfes = []
    for window_start_position in window_start_positions:
      mfes += [RNA.fold(seq[window_start_position:(window_start_position +
                                                   window_size)])[1]]
    # Return the min, median, max, and mean
    return min(mfes), np.median(mfes), max(mfes), sum(mfes)/len(mfes)
      
def get_mfe(df, seq_column="exonic_3UTR_seq", num_workers=8,
            window_size=100, step_size=50):
  ddata = dd.from_pandas(df, npartitions=num_workers)
  mfe_stats = (ddata
               .map_partitions(lambda df: 
                               df.apply(lambda row:
                                        mfe_sliding_window(row[seq_column],
                                                           window_size,
                                                           step_size),
                                        axis=1
                                       )
                               )
                 .compute(scheduler='processes')
                 .apply(pd.Series, index=['window_mfe_min', 'window_mfe_median', 
                                          'window_mfe_max', 'window_mfe_mean']))

  return pd.concat([df, mfe_stats], axis=1)
