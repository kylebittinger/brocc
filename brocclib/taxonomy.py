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
        "superkingdom", "kingdom", "phylum", "class",
        "order", "family", "genus", "species"
    ]
    standard_rank_idx = dict(
        (rank, idx) for idx, rank in enumerate(ranks))

    def __init__(self, taxa):
        self.indexed_lineage = [
            (name, self.standard_rank_idx.get(rank)) for name, rank in taxa]

    def get_taxon(self, rank):
        rank_idx = self.ranks.index(rank)
        for name, idx in self.indexed_lineage:
            if idx is None:
                continue
            elif idx == rank_idx:
                return name
            elif idx > rank_idx:
                return "{0} ({1})".format(name, rank)
        return None

    def get_standard_taxa(self, rank):
        rank_idx = self.ranks.index(rank)
        slice_idx = rank_idx + 1
        ranks_up_to = self.ranks[:slice_idx]
        for r in ranks_up_to:
            yield self.get_taxon(r)

    def get_all_taxa(self, rank):
        rank_idx = self.ranks.index(rank)
        for name, idx in self.indexed_lineage:
            if idx is None:
                yield name
            elif idx < rank_idx:
                yield name
            if idx == rank_idx:
                yield name
                break
            if idx > rank_idx:
                # If we've skipped over the rank of interest, yield a
                # placeholder and stop
                yield self.get_taxon(rank)
                break

    @staticmethod
    def is_generic_name(name):
        return (name == "environmental samples") or \
            name.startswith("unclassified")

    def is_generic(self, rank):
        rank_idx = self.ranks.index(rank)
        last_idx = 0
        for name, idx in self.indexed_lineage:
            if idx is None:
                if self.is_generic_name(name):
                    if last_idx < rank_idx:
                        return True
            else:
                if idx >= rank_idx:
                    return False
                last_idx = idx
        return False
