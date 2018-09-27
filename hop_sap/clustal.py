'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import os

from synbiochem.utils import seq_utils

import numpy as np
import pandas as pd


def align(df, out_dir):
    '''Align.'''
    result_file = os.path.join(out_dir, 'clustal.fasta')

    # Perform Clustal Omega alignment:
    seqs = pd.Series({k: v
                      for k, v in dict(zip(df['Entry name'],
                                           df.Sequence)).iteritems()
                      if not isinstance(v, np.float)})

    seqs_df = \
        pd.DataFrame.from_dict(seq_utils.do_clustal(seqs,
                                                    result_file=result_file),
                               columns=['align_seq'],
                               orient='index')

    # seqs_df.index = seqs_df.index.map(int)

    return df.merge(seqs_df, how='left',
                    left_on='Entry name', right_index=True)
