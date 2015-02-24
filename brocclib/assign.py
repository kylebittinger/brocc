from __future__ import division

from brocclib.Taxon import Lineage, NoLineage

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


class Assignment(object):
    def __init__(self, query_id, winning_candidate, total_votes, num_generic):
        self.id = query_id
        self.winning_candidate = winning_candidate
        self.total_votes = total_votes
        self.num_generic = num_generic

    def format_for_full_taxonomy(self):
        lineage = ';'.join(self.winning_candidate.all_taxa)
        return "%s\t%s\n" % (self.id, lineage)

    def format_for_standard_taxonomy(self):
        lineage = ';'.join(self.winning_candidate.standard_taxa)
        return "%s\t%s\n" % (self.id, lineage)

    def format_for_log(self):
        lineage = ';'.join(self.winning_candidate.all_taxa)
        return "%s\t%s\t%s\t%s\t%s\t%s\n" % (
            self.id, self.winning_candidate.votes, self.total_votes,
            self.num_generic, self.winning_candidate.rank, lineage)


class NoAssignment(object):
    """Null object representing no assignment, with message."""
    def __init__(self, name, message):
        self.name = name
        self.message = message

    def vote(self, *args, **kwargs):
        pass

    @property
    def winning_candidate(self):
        return self.message

    @property
    def id(self):
        return self.name

    def format_for_full_taxonomy(self):
        return "%s\t%s\n" % (self.id, self.message)

    format_for_standard_taxonomy = format_for_full_taxonomy
    format_for_log = format_for_full_taxonomy


class Assigner(object):
    ranks = [
        "species", "genus", "family", "order",
        "class", "phylum", "kingdom", "domain",
        ]

    def __init__(self, min_cover, species_min_id, genus_min_id, min_id,
                 consensus_thresholds, max_generic, taxa_db):
        self.consensus_thresholds = consensus_thresholds
        self.max_generic = max_generic
        self.min_id = min_id
        self.min_cover = min_cover
        # Minimum identity for consideration at each rank
        self.rank_min_ids = [
            species_min_id, genus_min_id, min_id, min_id,
            min_id, min_id, min_id, min_id,
            ]
        self.taxa_db = taxa_db

    def _quality_filter(self, seq, hits):
        hits_to_keep = []
        num_low_coverage = 0
        for hit in hits:
            identity_is_ok = (hit.pct_id >= self.min_id)
            coverage_is_ok = (hit.coverage(seq) >= self.min_cover)

            if identity_is_ok and coverage_is_ok:
                hits_to_keep.append(hit)

            elif identity_is_ok and not coverage_is_ok:
                num_low_coverage += 1

        frac_low_coverage = num_low_coverage / len(hits)
        return hits_to_keep, frac_low_coverage

    def assign(self, name, seq, hits):
        if not hits:
            return NoAssignment(name, "No hits found in database")

        hits_to_keep, frac_low_coverage = self._quality_filter(seq, hits)

        if frac_low_coverage > .9:
            return NoAssignment(name, "Abundance of low coverage hits: possible chimera")

        if not hits_to_keep:
            return NoAssignment(name, "All BLAST hits were filtered for low quality.")

        return self.vote(name, seq, hits_to_keep)

    def _retrieve_lineage(self, hit):
        taxid = self.taxa_db.get_taxon_id(hit.gi)
        if taxid is None:
            return NoLineage()

        raw_lineage = self.taxa_db.get_lineage(taxid)
        if raw_lineage is None:
            return NoLineage()

        return Lineage(raw_lineage)

    def vote(self, name, seq, hits):
        # Sort hits by percent ID.  This affects the way that ties are broken.
        hits.sort(reverse=True, key=lambda x: x.pct_id)
        hits_lineage = [(hit, self._retrieve_lineage(hit)) for hit in hits]
        for rank in self.ranks:
            a = self.vote_at_level(name, rank, hits_lineage)
            if a is not None:
                return a
        return NoAssignment(
            name, "Could not find consensus at domain level. No classification.")

    def vote_at_level(self, query_id, rank, db_hits):
        '''Votes at a given level (rank) of the taxonomy.'''
        rank_idx = self.ranks.index(rank)
        min_pct_id = self.rank_min_ids[rank_idx]
        consensus_threshold = self.consensus_thresholds[rank_idx]

        # Cast votes and count generic taxa
        candidates = dict()
        num_generic = 0
        for hit, lineage in db_hits:
            if hit.pct_id <= min_pct_id:
                continue
            taxon = lineage.get_taxon(rank)
            if taxon is None:
                continue
            if taxon not in candidates:
                candidates[taxon] = Candidate(lineage, rank)
            candidates[taxon].votes += 1

            if lineage.classified is False:
                num_generic += 1

        # Immediately return if there are no candidates
        if len(candidates) == 0:
            return None

        total_votes = sum(c.votes for c in candidates.values())

        # Do not count the votes for generic candidates in the total,
        # when the proportion of generic votes is allowable.  There
        # are some issues here with generic taxa, though the software
        # normally does the right thing.
        if (num_generic / total_votes) < self.max_generic:
            total_votes = total_votes - num_generic

        # Sort candidates by total votes and pick the winner.
        sorted_candidates = candidates.items()
        sorted_candidates.sort(reverse=True, key=lambda (taxon, c): c.votes)
        (name, winner) = sorted_candidates[0]

        # If there is a 50-50 tie for first place and the winner is a
        # bad classification, consider the runner up.
        if (winner.legit_taxa is False) and (len(sorted_candidates) > 1):
            (name, winner) = sorted_candidates[1]
            # If both top hits are low quality, there is no winner
            if winner.legit_taxa is False:
                return None

        # Only one low quality hit, no winner.
        elif winner.legit_taxa is False and len(candidates) == 1:
            return None

        if (winner.votes / total_votes) > consensus_threshold:
            return Assignment(query_id, winner, total_votes, num_generic)
        else:
            return None
