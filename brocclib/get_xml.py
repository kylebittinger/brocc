import json
import urllib2, os, StringIO
from xml.etree import ElementTree as ET
import logging


class NcbiEutils(object):
    def __init__(self, cache_fp=None):
        self.cache_fp = cache_fp
        self.lineages = {}
        self.taxon_ids = {}
        self._fresh = True

    def get_lineage(self, taxon_id):
        if taxon_id not in self.lineages:
            self._fresh = False
            self.lineages[taxon_id] = get_lineage(taxon_id)
        return self.lineages[taxon_id]

    def get_taxon_id(self, acc):
        if acc not in self.taxon_ids:
            self._fresh = False
            self.taxon_ids[acc] = get_taxid(acc)
        return self.taxon_ids[acc]

    def load_cache(self):
        # Do nothing if there is no cache file
        if self.cache_fp is None:
            return None
        if os.path.exists(self.cache_fp):
            with open(self.cache_fp) as f:
                data = json.load(f)
                self.lineages = dict(
                    (x, dict(y)) for x, y in data["lineages"])
                self.taxon_ids = dict(data["taxon_ids"])

    def save_cache(self):
        # Do nothing if there is no cache file
        if self.cache_fp is None:
            return None
        # Do nothing if no lookups have been performed.
        if self._fresh:
            return None

        lineages = list(
            (x, list(sorted(y.items()))) for x, y in self.lineages.items())
        lineages.sort()

        taxon_ids = list(self.taxon_ids.items())
        taxon_ids.sort()

        data = {
            "lineages": lineages,
            "taxon_ids": taxon_ids,
            }
        with open(self.cache_fp, "w") as f:
            json.dump(data, f, indent=2, separators=(',', ': '))


def get_taxon_from_xml(xml_string):
    taxon_dict = {}
    xml_string_new = "".join(
        [s for s in xml_string.splitlines(True) if s.strip("\r\n")])
    tree = ET.XML(xml_string_new)
    lineage_elem = tree.find('Taxon/LineageEx')
    if lineage_elem is None:
        raise ValueError("No lineage info found in XML:\n" + xml_string)
    for elem in list(lineage_elem):
        rank = elem.find('Rank').text
        name = elem.find('ScientificName').text
        taxon_dict[rank] = name

    # Include lowest rank in lineage
    rank = tree.find('Taxon/Rank').text
    if rank not in taxon_dict:
        taxon_dict[rank] = tree.find('Taxon/ScientificName').text

    # Also include the lineage as a string
    taxon_dict['Lineage'] = tree.find('Taxon/Lineage').text

    return taxon_dict


def _get_xml_from_html(html_response):
    if not html_response.strip():
        raise ValueError("Empty HTML response")
    html_tree = ET.fromstring(html_response)
    html_body = html_tree.find('body')
    if html_body is None:
        # Response is already XML
        return html_response
    else:
        # Need to find XML string embedded in HTML
        s = html_body.find("pre").text
        s = s.replace("&lt;", "<")
        s = s.replace("&gt;", ">")
        # this has to be last:
        s = s.replace("&amp;", "and")
        return s


def get_lineage(taxid):
    num_tries = 0 #numter of times db connection was attempted
    while num_tries < 5:
        try:        #watch out for db connection time out
            url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=' + taxid + '&rettype=xml'
            xml = urllib2.urlopen(url)
            xml_str = xml.read()
            xml_string_from_html = _get_xml_from_html(xml_str)
            taxon_dict = get_taxon_from_xml(xml_string_from_html)
            return taxon_dict
        except urllib2.HTTPError as e:
            logging.info(str(e))
            if e.code == 400:
                return None
            else:
                num_tries += 1
        except urllib2.URLError as e:
            num_tries += 1
            print e
            print "Database connectiom timed out for taxon", taxid, "Will retry"
        except Exception as e:
            num_tries += 1
            print e, "Will retry"
    if num_tries == 5:
        print "could not connect to db for taxon" + str(taxid)
        print str(taxid) + " will not be considered"
        return None


def url_open(url, max_tries=5):
    for n in xrange(max_tries):
        try:
            return urllib2.urlopen(url)
        except urllib2.HTTPError as e:
            # Don't keep trying if you gave a bad request
            if e.code == 400:
                raise e
            logging.debug("Retrying URL %s (attempt %s)" % (url, n))
        except urllib2.URLError:
            logging.debug("Retrying URL %s (attempt %s)" % (url, n))
    raise urllib2.URLError("Could not open URL %s (%s attempts)" % (url, n))


def get_taxid(acc):
    url_fmt = (
        "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?"
        "dbfrom=nucleotide&db=taxonomy&id={}")
    url = url_fmt.format(acc)
    xpath = ".//Link/Id"
    try:
        response = url_open(url)
        xml = ET.parse(response)
        elem = xml.find(xpath)
        if elem is None:
            logging.debug("Accession %s: Xpath %s not found in XML" % (acc, xpath))
            return None
        return elem.text
    except Exception as e:
        logging.info("Accession %s: %s" % (acc, e))
        return None
