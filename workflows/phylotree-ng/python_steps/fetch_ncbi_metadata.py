import argparse
import json
import logging
import requests
import time
import xml.etree.ElementTree as ET
from urllib.parse import urlencode
from typing import Iterable, Optional, TypedDict


class Metadata(TypedDict, total=False):
    name: Optional[str]
    country: Optional[str]
    collection_date: Optional[str]


def get_accession_metadata(accession_id: str):
    '''
    Retrieve metadata of an NCBI accession (e.g. name, country, collection date)
    TODO: Put this data in S3 instead and get it from there.
    '''
    accession_metadata: Metadata = {}
    try:
        baseurl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        output = requests.get(f"{baseurl}/esearch.fcgi?" + urlencode({
            "db": "nuccore",
            "term": accession_id,
            "usehistory": "y"
        })).text
        # sleep to avoid throttling
        time.sleep(1)

        root = ET.fromstring(output)
        web = root.find('WebEnv').text
        key = root.find('QueryKey').text
        genbank_xml = requests.get(f"{baseurl}/efetch.fcgi?" + urlencode({
            "db": "nuccore",
            "query_key": key,
            "WebEnv": web,
            "rettype": "gb",
            "retmode": "xml"
        })).text
        # sleep to avoid throttling
        time.sleep(1)

        root = ET.fromstring(genbank_xml).find('GBSeq')
        if not root:
            logging.warn(f"Failed to fetch metadata for accession ID: {accession_id}")
            return accession_metadata
        accession_metadata['name'] = root.find('GBSeq_definition').text
        qualifiers_needed = {'country', 'collection_date'}
        for entry in root.find('GBSeq_feature-table')[0].find('GBFeature_quals'):
            if all(key in accession_metadata for key in qualifiers_needed):
                break
            for key in qualifiers_needed - accession_metadata.keys():
                if entry.find('GBQualifier_name').text == key:
                    accession_metadata[key] = entry.find('GBQualifier_value').text
        return accession_metadata
    except Exception as e:
        logging.warn(f"Fetching metadata for accession ID {accession_id} failed with error: {e}", exec_info=True)
        return accession_metadata


def main(reference_accession_ids: Iterable[str], output_ncbi_metadata: str):
    with open(output_ncbi_metadata, "w") as f:
        json.dump({a: get_accession_metadata(a) for a in reference_accession_ids}, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-accession-ids")
    parser.add_argument("--output-ncbi-metadata")
    args = parser.parse_args()

    with open(args.reference_accession_ids) as f:
        reference_accession_ids: Iterable[str] = json.load(f)

    main(reference_accession_ids, args.output_ncbi_metadata)
