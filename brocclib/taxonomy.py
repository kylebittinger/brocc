'''
Created on Aug 29, 2011
@author: Serena
'''

class NoLineage(object):
    def get_taxon(self, rank):
        return None


class Lineage(object):
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
