'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
from synbiochem.utils import seq_utils

import numpy as np
import pandas as pd


def align(df):
    '''Align.'''

    # Perform Clustal Omega alignment:
    seqs = pd.Series({k: v
                      for k, v in df.to_dict()['Sequence'].iteritems()
                      if not isinstance(v, np.float)})

    seqs_df = pd.DataFrame.from_dict(seq_utils.do_clustal(seqs),
                                     columns=['align_seq'],
                                     orient='index')

    seqs_df.index = seqs_df.index.map(int)

    return df.join(seqs_df)
