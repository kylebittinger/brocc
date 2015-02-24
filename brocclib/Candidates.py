'''
Created on Aug 29, 2011
@author: Serena, Kyle
'''

class Candidate(object):
    """Candidate assignment"""
    def __init__(self, lineage, rank):
        self.votes = 0
        self.lineage = lineage
        self.rank = rank
        
        is_high_rank = rank in ["phylum", "kingdom", "domain"]
        is_descended_from_missing_taxon = any("(" in h for h in self.all_taxa)

        if is_high_rank and not is_descended_from_missing_taxon:
            self.legit_taxa = True
        else:
            self.legit_taxa = lineage.classified
        
    @property
    def all_taxa(self):
        return self.lineage.get_all_taxa(self.rank)

    @property
    def standard_taxa(self):
        return self.lineage.get_standard_taxa(self.rank)
