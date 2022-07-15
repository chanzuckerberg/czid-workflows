import argparse
import datetime
import logging
import multiprocessing
import os
import re
import time
from functools import update_wrapper

import marisa_trie
import pandas as pd


logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')


def decorator(d):
    "Make function d a decorator: d wraps a function fn."
    def _d(fn):
        return update_wrapper(d(fn), fn)
    update_wrapper(_d, d)
    return _d


@decorator
def timeit(f):
    """time a function, used as decorator"""
    def new_f(*args, **kwargs):
        bt = time.time()
        r = f(*args, **kwargs)
        et = time.time()
        delta_t = datetime.timedelta(seconds=(et - bt))
        logging.info("Time spent on {0}: {1}".format(f.__name__, delta_t))
        return r
    return new_f


def backup_file(f):
    """
    Back up a file, old_file will be renamed to #old_file.n#, where n is a
    number incremented each time a backup takes place
    """
    if os.path.exists(f):
        dirname = os.path.dirname(f)
        basename = os.path.basename(f)
        count = 1
        rn_to = os.path.join(
            dirname, '#' + basename + '.{0}#'.format(count))
        while os.path.exists(rn_to):
            count += 1
            rn_to = os.path.join(
                dirname, '#' + basename + '.{0}#'.format(count))
        logging.info("Backing up {0} to {1}".format(f, rn_to))
        os.rename(f, rn_to)
        return rn_to
    else:
        logging.warning('{0} doesn\'t exist'.format(f))


def parse_args():
    parser = argparse.ArgumentParser(
        description=('This script converts NCBI taxonomy (taxdump) into '
                     'lineages, save the information in csv.gz format'))

    parser.add_argument(
        '--nodes-file', required=True,
        help='NCBI taxonomy path/to/taxdump/nodes.dmp')

    parser.add_argument(
        '--names-file', required=True,
        help='NCBI taxonomy path/to/taxdump/names.dmp')

    parser.add_argument(
        '-o', '--output-prefix', default='ncbi_lineages',
        help='will output lineage name information in [output_prefix].csv.gz')

    parser.add_argument(
        '--names-output-prefix', default='ncbi_names',
        help='will output scientific-name information in [names_output_prefix].csv.gz')

    parser.add_argument(
        '--taxid-lineages-output-prefix', default='ncbi_taxid_lineages',
        help='will output lineage taxon-ID information in [taxid_lineages_output_prefix].csv.gz')

    parser.add_argument(
        '--name-lineages-output-prefix', default='ncbi_name_lineages',
        help='will output lineage name information in [name_lineages_output_prefix].csv.gz')

    args = parser.parse_args()
    return args


def strip(str_):
    '''
    :param str_: a string
    '''
    return str_.strip()


@timeit
def load_nodes(nodes_file):
    '''
    load nodes.dmp and convert it into a pandas.DataFrame
    '''
    df = pd.read_csv(nodes_file, sep='|', header=None, index_col=False,
                     names=[
                         'tax_id',
                         'parent_tax_id',
                         'rank',
                         'embl_code',
                         'division_id',
                         'inherited_div_flag',
                         'genetic_code_id',
                         'inherited_GC__flag',
                         'mitochondrial_genetic_code_id',
                         'inherited_MGC_flag',
                         'GenBank_hidden_flag',
                         'hidden_subtree_root_flag',
                         'comments'
                     ])

    # To get rid of flanking tab characters
    df['rank'] = df['rank'].apply(strip)
    df['embl_code'] = df['embl_code'].apply(strip)
    df['comments'] = df['comments'].apply(strip)
    return df


@timeit
def load_names(names_file, name_class='scientific name'):
    '''
    load names.dmp and convert it into a pandas.DataFrame
    '''
    df = pd.read_csv(names_file, sep='|', header=None, index_col=False,
                     names=[
                         'tax_id',
                         'name_txt',
                         'unique_name',
                         'name_class'
                     ])
    df['name_txt'] = df['name_txt'].apply(strip)
    df['unique_name'] = df['unique_name'].apply(strip)
    df['name_class'] = df['name_class'].apply(strip)

    sci_df = df[df['name_class'] == name_class]
    sci_df.reset_index(drop=True, inplace=True)
    return sci_df


def to_name_dict(lineage):
    """
    convert the lineage from a list of tuples in the form of
    [
        (tax_id1, rank1, name_txt1),
        (tax_id2, rank2, name_txt2),
        ...
    ]
    to a dictionary of taxon names
    """
    dd = {}
    num_re = re.compile('[0-9]+')
    len_lineage = len(lineage)
    for k, __ in enumerate(lineage):
        tax_id, rank, name_txt = __
        # use the last rank as the tax_id, whatever it is, genus or species.
        if k == len_lineage - 1:
            dd['tax_id'] = tax_id

        # e.g. there could be multiple 'no rank'
        numbered_rank = rank
        while numbered_rank in dd:
            search = num_re.search(numbered_rank)
            if search is None:
                count = 1
            else:
                count = int(search.group()) + 1
            numbered_rank = '{0}{1}'.format(rank, count)
        dd[numbered_rank] = name_txt
    return dd


def to_taxid_dict(lineage):
    """
    convert the lineage from a list of tuples in the form of
    [
        (tax_id1, rank1, name_txt1),
        (tax_id2, rank2, name_txt2),
        ...
    ]
    to a dictionary of taxids
    """
    dd = {}
    num_re = re.compile('[0-9]+')
    len_lineage = len(lineage)
    for k, __ in enumerate(lineage):
        tax_id, rank, name_txt = __
        # use the last rank as the tax_id, whatever it is, genus or species.
        if k == len_lineage - 1:
            dd['tax_id'] = int(tax_id)

        # e.g. there could be multiple 'no rank'
        numbered_rank = rank
        while numbered_rank in dd:
            search = num_re.search(numbered_rank)
            if search is None:
                count = 1
            else:
                count = int(search.group()) + 1
            numbered_rank = '{0}{1}'.format(rank, count)
        dd[numbered_rank] = int(tax_id)
    return dd


def find_lineage(tax_id):
    if tax_id % 50000 == 0:
        logging.debug('working on tax_id: {0}'.format(tax_id))
    lineage = []
    while True:
        rec = TAXONOMY_DICT[tax_id]
        lineage.append((rec['tax_id'], rec['rank'], rec['name_txt']))
        tax_id = rec['parent_tax_id']

        if tax_id == 1:
            break

    # reverse results in lineage of Kingdom => species, this is helpful for
    # to_dict when there are multiple "no rank"s
    lineage.reverse()
    return to_name_dict(lineage), to_taxid_dict(lineage)


def process_lineage_dd(lineage_dd):
    dd_for_df = dict(zip(range(len(lineage_dd)), lineage_dd))
    lineages_df = pd.DataFrame.from_dict(dd_for_df, orient='index')
    return lineages_df.sort_values('tax_id')


def write_output(output_prefix, output_name_log, df, cols=None, undef_taxids=None):
    output = os.path.join('{0}.csv.gz'.format(output_prefix))
    logging.info("writing %s to %s" % (output_name_log, output))
    if undef_taxids and cols:
        for col in undef_taxids.keys():
            df[[col]] = df[[col]].fillna(value=undef_taxids[col])
        # filling remaing na values as 0
        df[cols] = df[cols].fillna(0)
        df[cols] = df[cols].astype(int)
    df.to_csv(output, index=False, columns=cols, compression='gzip')


def generate_name_output(nodes_df, names_file, name_class):
    names_df = load_names(names_file, name_class)
    df = nodes_df.merge(names_df, on='tax_id')
    df = df[['tax_id', 'parent_tax_id', 'rank', 'name_txt']]
    df.reset_index(drop=True, inplace=True)
    logging.info('# of tax ids: {0}'.format(df.shape[0]))
    df.info()
    return df


@timeit
def generate_lineage_outputs(df, taxid_lineages_output_prefix, name_lineages_output_prefix):
    def lineage_by_taxid():
        # example item: (16, {'parent_tax_id': 32011, 'name_txt': 'Methylophilus', 'rank': 'genus', 'tax_id': 16})
        global TAXONOMY_DICT
        logging.info('generating TAXONOMY_DICT...')
        TAXONOMY_DICT = dict(zip(df.tax_id.values, df.to_dict('records')))

        ncpus = multiprocessing.cpu_count()
        logging.info('found {0} cpus, and will use all of them to find lineages '
                     'for all tax ids'.format(ncpus))
        pool = multiprocessing.Pool(ncpus)
        name_lineages_dd, taxid_lineages_dd = zip(*pool.map(find_lineage, df.tax_id.values))  # take about 18G memory
        pool.close()

        logging.info('generating lineage-by-name output...')
        name_lineages_df = process_lineage_dd(name_lineages_dd)
        name_lineages_df.columns = name_lineages_df.columns.str.replace('no rank', 'no_rank')
        write_output(name_lineages_output_prefix, "name lineages", name_lineages_df,
                     ['tax_id', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'] +
                     [col for col in name_lineages_df if col.startswith('no_rank')])

        logging.info('generating lineage-by-taxid output...')
        taxid_lineages_df = process_lineage_dd(taxid_lineages_dd)
        taxid_lineages_df.columns = taxid_lineages_df.columns.str.replace(' ', '_')
        undef_taxids = {'species': -100,
                        'genus': -200,
                        'family': -300,
                        'order': -400,
                        'class': -500,
                        'phylum': -600,
                        'kingdom': -650,
                        'superkingdom': -700}
        write_output(taxid_lineages_output_prefix, "taxid lineages", taxid_lineages_df,
                     ['tax_id', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'] +
                     [col for col in taxid_lineages_df if col.startswith('no_rank')],
                     undef_taxids=undef_taxids)

        logging.info('writing lineage-by-taxid ...')
        for _, row in taxid_lineages_df.iterrows():
            yield (str(int(row['tax_id'])), (int(row['species']), int(row['genus']), int(row['family'])))
    marisa_trie.RecordTrie("III", lineage_by_taxid()).save()


def main():
    args = parse_args()

    logging.info('PART I: name outputs')
    nodes_df = load_nodes(args.nodes_file)
    logging.info(' * get scientific names')
    scientific_df = generate_name_output(nodes_df, args.names_file, 'scientific name')
    logging.info(' * get common names')
    common_df = generate_name_output(nodes_df, args.names_file, 'genbank common name')
    logging.info(' * merge names')
    df = scientific_df.merge(common_df, 'left', on='tax_id', suffixes=('', '_common'))
    write_output(args.names_output_prefix, "names", df, ['tax_id', 'name_txt', 'name_txt_common'])

    logging.info('PART II: lineage and scientific name outputs')
    generate_lineage_outputs(scientific_df, args.taxid_lineages_output_prefix, args.name_lineages_output_prefix)


if __name__ == "__main__":
    main()
