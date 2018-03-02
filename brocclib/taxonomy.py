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
        "class", "phylum", "kingdom", "superkingdom",
        ]

    def __init__(self, dictionary):
        self.store = dictionary
        self.full_lineage = self.store["LineageWithRanks"]
        self.standard_taxa = self._create_standard_taxa(self.full_lineage)

    @classmethod
    def _create_standard_taxa(cls, lineage):
        std_taxa = dict(
            (rank, name) for name, rank in lineage if rank in cls.ranks)
        fill_in_val = None
        for rank in cls.ranks:
            if rank in std_taxa:
                fill_in_val = std_taxa[rank]
            else:
                std_taxa[rank] = "{0} ({1})".format(fill_in_val, rank)
        return std_taxa

    def get_standard_taxa(self, rank):
        for r in reversed(self.ranks):
            t = self.get_taxon(r)
            if t is not None:
                yield t
            if rank == r:
                break

    def get_all_taxa(self, rank):
        for taxon_name, taxon_rank in self.full_lineage:
            yield taxon_name
            if taxon_rank == rank:
                break

    def get_taxon(self, rank):
        return self.standard_taxa.get(rank)

    def is_generic(self, rank):
        return is_generic(self.get_taxon(rank))
