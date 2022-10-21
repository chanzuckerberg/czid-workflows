import json
import threading
import re
import time
import wikipedia
from Bio import Entrez
import idseq_dag.util.log as log
import idseq_dag.util.s3 as s3
from idseq_dag.engine.pipeline_step import PipelineStep
from idseq_dag.util.trace_lock import TraceLock


# This does not run from the actual pipeline.  It's an index generation step.
class PipelineStepFetchTaxInfo(PipelineStep):
    '''
        fetch tax info based on a list
    '''
    def run(self):
        '''
            1. fetch the taxid -> wikipedia link mapping
            2. fetch wikipedia content
            3. store everything
        '''
        taxid_list = self.input_files_local[0][0]
        (taxid2wiki, taxid2desc) = self.output_files_local()

        taxid2wikidict = {}
        Entrez.email = self.additional_attributes.get("entrez_email", "idseq-tech@chanzuckerberg.com")
        num_threads = self.additional_attributes.get("threads", 16)
        batch_size = self.additional_attributes.get("batch_size", 100)
        namecsv = self.additional_files.get("taxon2name")
        id2namedict = {}
        if namecsv:
            # This is fetching a reference without fetch_reference;  but ok because does not run from the actual pipeline
            namecsvf = s3.fetch_from_s3(namecsv, "/mnt/idseq/ref")
            with open(namecsvf, 'r') as namef:
                for line in namef:
                    fields = line.rstrip().split(",")
                    id2namedict[fields[0]] = fields[1]

        # This is fetching a reference without fetch_reference and doing a presence check;  but ok because does not run from the actual pipeline
        if s3.check_s3_presence(self.s3_path(taxid2wiki)):
            # generated
            taxid2wiki = s3.fetch_from_s3(self.s3_path(taxid2wiki), taxid2wiki)
            with open(taxid2wiki, "r") as taf:
                for line in taf:
                    (key, val) = line.rstrip("\n").split("\t")
                    taxid2wikidict[key] = val
        else:
            self.fetch_ncbi_wiki_map(num_threads, batch_size, taxid_list, taxid2wikidict)

        # output dummay for actual wiki content for now
        taxid2wikicontent = {}
        self.fetch_wiki_content(num_threads * 4, taxid2wikidict,
                                taxid2wikicontent, id2namedict)

        with open(taxid2desc, 'w') as desc_outputf:
            json.dump(taxid2wikicontent, desc_outputf)

        # output the taxid 2 wikiurl data
        with open(taxid2wiki, 'w') as taxidoutf:
            for taxid, wikiurl in taxid2wikidict.items():
                if wikiurl == "":
                    pageid = taxid2wikicontent.get(taxid, {}).get('pageid', None)
                    if pageid:
                        wikiurl = f"http://en.wikipedia.org/wiki/index.html?curid={pageid}"
                taxidoutf.write(f"{taxid}\t{wikiurl}\n")

    @staticmethod
    def fetch_wiki_content(num_threads, taxid2wikidict, taxid2wikicontent, id2namedict):
        ''' Fetch wikipedia content based on taxid2wikidict '''
        threads = []
        semaphore = threading.Semaphore(num_threads)
        mutex = TraceLock("fetch_wiki_content", threading.RLock())
        for taxid, url in taxid2wikidict.items():
            m = re.search(r"curid=(\d+)", url)
            pageid = None
            if m:
                pageid = m[1]
            name = id2namedict.get(taxid)
            if pageid or name:
                semaphore.acquire()
                t = threading.Thread(
                    target=PipelineStepFetchTaxInfo.
                    get_wiki_content_for_page,
                    args=[taxid, pageid, name, taxid2wikicontent, mutex, semaphore]
                )
                t.start()
                threads.append(t)
        for t in threads:
            t.join()

    @staticmethod
    def get_wiki_content_for_page(taxid, pageid, taxname, taxid2wikicontent, mutex, semaphore, max_attempt=3):
        ''' Fetch wiki content for pageid '''
        for attempt in range(max_attempt):
            try:
                page = None
                if pageid:
                    log.write(f"fetching wiki {pageid} for {taxid}")
                    page = wikipedia.page(pageid=pageid)
                elif taxname:
                    search_results = wikipedia.search(taxname)
                    if len(search_results) > 0:
                        wikiname = str(search_results[0])
                        if taxname.lower() == wikiname.lower():
                            page = wikipedia.page(wikiname)
                    if not page:
                        # query the page directly
                        try:
                            page = wikipedia.page(taxname.replace(" ", "_"))
                        except:
                            page = None

                if page:
                    output = {
                        "pageid": page.pageid,
                        "description": page.content[:1000],
                        "title": page.title,
                        "summary": page.summary
                    }
                    with mutex:
                        taxid2wikicontent[taxid] = output
                break
            except:
                log.write(f"having trouble fetching {taxid} wiki {pageid} attempt {attempt}")
        semaphore.release()

    @staticmethod
    def fetch_ncbi_wiki_map(num_threads, batch_size, taxid_list, taxid2wikidict):
        ''' Use Entrez API to fetch taxonid -> wikipedia page mapping '''
        threads = []
        semaphore = threading.Semaphore(num_threads)
        mutex = TraceLock("fetch_ncbi_wiki_map", threading.RLock())
        batch = []
        with open(taxid_list, 'r') as taxf:
            for line in taxf:
                taxid = line.rstrip()
                if taxid == 'taxid':
                    continue  # header
                batch.append(taxid)
                if len(batch) >= batch_size:
                    semaphore.acquire()
                    t = threading.Thread(
                        target=PipelineStepFetchTaxInfo.
                        get_taxid_mapping_for_batch,
                        args=[batch, taxid2wikidict, mutex, semaphore]
                    )
                    t.start()
                    threads.append(t)
                    batch = []
        if len(batch) > 0:
            semaphore.acquire()
            t = threading.Thread(
                target=PipelineStepFetchTaxInfo.
                get_taxid_mapping_for_batch,
                args=[batch, taxid2wikidict, mutex, semaphore]
            )
            t.start()
            threads.append(t)
        for t in threads:
            t.join()

    @staticmethod
    def get_taxid_mapping_for_batch(taxids, taxid2wikidict, mutex, semaphore, max_attempt=3):
        ''' Get wiki mapping for a list of taxids '''
        taxid_str = ",".join(taxids)
        log.write(f"fetching batch {taxid_str}")
        for attempt in range(max_attempt):
            try:
                handle = Entrez.elink(dbfrom="taxonomy", id=taxid_str, cmd="llinks")
                record = Entrez.read(handle)
                handle.close()

                parsed = {}
                results = record[0]['IdUrlList']['IdUrlSet']
                for result in results:
                    taxid = result['Id']
                    wikiurl = ""
                    for link in result['ObjUrl']:
                        url = str(link['Url'])
                        if re.search('wikipedia.org', url):
                            wikiurl = url
                            break
                    parsed[taxid] = wikiurl
                break
            except:
                log.write(f"failed batch attempt {attempt}")
                time.sleep(5)
        semaphore.release()
        with mutex:
            taxid2wikidict.update(parsed)

    def count_reads(self):
        pass
