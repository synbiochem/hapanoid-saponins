'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=wrong-import-order
import csv
import os
import sys
from urllib import urlretrieve

from hap_sap import clustal
import pandas as pd


def get_seqs(in_filename, out_dir):
    '''Get Uniprot sequences.'''
    df = _get_data(in_filename)
    uniprot_df = _get_uniprot_data(df, out_dir)

    return df.merge(uniprot_df,
                    how='left',
                    left_on='Genbank', right_on='Cross-reference (EMBL)')


def _get_data(filename):
    '''Get data.'''
    return pd.read_csv(filename, encoding='latin1')


def _get_uniprot_data(df, out_dir):
    '''Get Uniprot data.'''
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    uniprot_csv = os.path.join(out_dir, 'uniprot.tsv')

    if not os.path.exists(uniprot_csv):
        genbank_ids = list(df['Genbank'].dropna())
        query = '+OR+'.join(genbank_ids)

        query_str = query + \
            '&format=tab' + \
            '&columns=id,entry name,protein names,organism,organism-id,' + \
            'ec,go,go-id,database(EMBL),sequence'

        url = 'http://www.uniprot.org/uniprot/?query=' + query_str

        urlretrieve(url, uniprot_csv)

    # Read Uniprot data into Dataframe:
    embl = 'Cross-reference (EMBL)'
    uniprot_df = pd.read_csv(uniprot_csv, sep='\t')
    return _split_data_frame_list(uniprot_df[uniprot_df[embl].notnull()],
                                  column=embl,
                                  sep=';')


def _split_data_frame_list(df, column, sep):
    '''Split list from column into individual rows.'''
    new_rows = []
    df.apply(_split_list_to_rows, axis=1, args=(new_rows, column, sep))
    return pd.DataFrame(new_rows)


def _split_list_to_rows(row, new_rows, column, sep):
    '''Split list into rows.'''
    for term in row[column].split(sep):
        new_row = row.to_dict()
        new_row[column] = term
        new_rows.append(new_row)


def main(args):
    '''main method.'''
    df = get_seqs(*args)
    df = clustal.align(df, args[1])
    df.to_csv(os.path.join(args[1], 'out.csv'),
              index=False, encoding='utf8', quoting=csv.QUOTE_NONNUMERIC)


if __name__ == '__main__':
    main(sys.argv[1:])
