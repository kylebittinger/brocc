'''
Created on Aug 29, 2011
@author: Serena
'''

GENERIC_TAXA = [
    "uncultured eukaryote",
    "ectomycorrhizal fungus",
    "endophytic basidiomycete",
    "fungal endophyte",
    "soil zygomycete",
    "uncultured Agaricomycetes",
    "Uncultured Agaricomycotina clone",
    "uncultured alveolate",
    "uncultured Ascomycota",
    "uncultured Basidiomycota",
    "uncultured ciliate",
    "uncultured compost fungus",
    "uncultured Dikarya",
    "uncultured Dothideomycetidae",
    "uncultured endophytic fungus",
    "uncultured eukaryote",
    "uncultured freshwater eukaryote",
    "uncultured fungal endophyte",
    "uncultured fungus",
    "uncultured Lecanoromycetida",
    "uncultured marine eukaryote",
    "uncultured marine picoeukaryote",
    "uncultured metazoan",
    "uncultured mycorrhizal fungus",
    "uncultured organism",
    "uncultured Pleosporales",
    "uncultured Pucciniomycotina",
    "uncultured root-associated fungus",
    "uncultured soil basidiomycete",
    "uncultured soil fungus",
    "uncultured stichotrichid",
    "uncultured zygomycete",
    "vouchered mycorrhizae",
    "Uncultured bacterium clone",
    "Uncultured organism clone",
    "uncultured Bacteria bacterium",
    "uncultivated bacterium",
    "unidentified bacterium",
    "unidentified bacteria",
    "unknown bacteria",
    "unclassified bacterium",
    "uncultured eubacterium",
    "Unknown eubacteria",
    "Unknown eubacterium",
    "unidentified eubacterium",
    "Uncultured Firmicutes bacterium clone",
    "Uncultured Bacteroidetes bacterium clone",
    ]

GENERIC_FLAGS = [
    "unclassified Fungi",
    "unclassified Basidiomycota",
    "unclassified Ascomycota",
    ]

GENERIC_WORDS = [
    # "environmental samples" also appears in lineage, but has no rank
    "uncultured", "unclassified", "unidentified", "fungal sp.",
    "fungal endophyte", "root associated fungus", "ectomycorrhizal fungus",
    "endophytic basidiomycete", "soil zygomycete", "vouchered", "fungal contaminant",
    "basidiomycete sp.", "Basidiomycota sp.", "Ascomycota sp.", "ascomycete strain",
]

# descended from "environmental samples" or "unclassified Fungi|Basidiomycota|XXXX"?
def is_generic(taxon_name):
    norm_taxon = taxon_name.lower()
    if any((gword in norm_taxon) for gword in GENERIC_WORDS):
        return True
    if taxon_name in GENERIC_TAXA:
        return True
    return False

class NoLineage(object):
    def get_taxon(self, rank):
        return None


class Lineage(object):
    generic_taxa = GENERIC_TAXA
    generic_flags = GENERIC_FLAGS
    ranks = [
        "species", "genus", "family", "order",
        "class", "phylum", "kingdom", "domain",
        ]

    def __init__(self, dictionary):
        self.store = dictionary

        self.species = self.store.get("species")

        self.classified = True
        if self.species in self.generic_taxa:
            self.classified = False
        if ("no rank" in self.store) and (self.store["no rank"] in self.generic_flags):
            self.classified = False

        ### FIXME: do not store taxa in attributes.
        self.genus = self.store.get("genus")
        if (self.genus is None) and (self.species is not None):
            self.genus = self.species + " (genus)"

        self.family = self.store.get("family")
        if (self.family is None) and (self.genus is not None):
            self.family = self.genus.split(" (")[0] + " (family)"

        self.order = self.store.get("order")
        if (self.order is None) and (self.family is not None):
            self.order = self.family.split(" (")[0] + " (order)"

        self.clas = self.store.get("class")
        if (self.clas is None) and (self.order is not None):
            self.clas = self.order.split(" (")[0] + " (class)"

        self.phylum = self.store.get("phylum")
        if (self.phylum is None) and (self.clas is not None):
            self.phylum = self.clas.split(" (")[0] + " (phylum)"

        self.kingdom = self.store.get("kingdom")
        if (self.kingdom is None) and (self.phylum is not None):
            self.kingdom = self.phylum.split(" (")[0] + " (kingdom)"

        self.domain = self.store.get("superkingdom")
        if (self.domain is None) and (self.kingdom is not None):
            self.domain = "Domain unknown for reference"
        ################## End FIXME

        self.full_lineage = self.store["Lineage"].split("; ")
        if self.species is not None:
            self.full_lineage.append(self.species)

    def get_standard_taxa(self, rank):
        for r in reversed(self.ranks):
            t = self.get_taxon(r)
            if t is not None:
                yield t
            if rank == r:
                break

    def get_all_taxa(self, rank):
        taxon = self.get_taxon(rank)
        for t in self.full_lineage:
            yield t
            if t == taxon:
                break

    def get_taxon(self, rank):
        if rank == "class":
            rank = "clas"
        return getattr(self, rank)
