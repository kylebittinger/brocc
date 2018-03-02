from __future__ import division
import collections
import logging

from brocclib.taxonomy import Lineage, NoLineage

'''
Created on Aug 29, 2011
@author: Serena, Kyle
'''

class AssignmentCandidate(object):
    def __init__(self, lineage, rank):
        self.votes = 0
        self.lineage = lineage
        self.rank = rank

    def is_placeholder(self):
        return "(" in self.lineage.get_taxon(self.rank)
        #is_high_rank = self.rank in ["phylum", "kingdom", "domain"]
        #is_descended_from_missing_taxon = any("(" in h for h in self.standard_taxa)
        #return not (is_high_rank and not is_descended_from_missing_taxon)

    @property
    def standard_taxa(self):
        return self.lineage.get_standard_taxa(self.rank)

    def to_string(self):
        return "{1} ({2} votes)".format(self.rank, self.lineage.get_taxon(self.rank), self.votes)


class Assignment(object):
    is_valid_assignment = True

    def __init__(self, query_id, winning_candidate=None,
                 rank=None, candidates=None, generics=None):
        self.query_id = query_id
        self.winning_candidate = winning_candidate
        self.rank = rank
        self.candidates = candidates or []
        self.generics = generics or []

    @property
    def total_votes(self):
        return sum(c.votes for c in self.candidates)

    @property
    def num_generic(self):
        return sum(n for n, taxon in self.generics)

    def format_for_standard_taxonomy(self):
        lineage = ';'.join(self.winning_candidate.standard_taxa)
        return "%s\t%s\n" % (self.query_id, lineage)

    def format_for_log(self):
        lineage = ';'.join(self.winning_candidate.standard_taxa)
        return "%s\t%s\t%s\t%s\t%s\t%s\n" % (
            self.query_id, self.winning_candidate.votes, self.total_votes,
            self.num_generic, self.winning_candidate.rank, lineage)

    def log_details(self):
        candidates2 = [c.to_string() for c in self.candidates if c.votes > 1]
        generics2 = ["{0} ({1} votes)".format(k, v) for v, k in self.generics]
        parts = [
            "",
            "** {0} voting for {1} **".format(self.query_id, self.rank),
            "Candidates (>1 vote): {0}".format(", ".join(candidates2)),
            "Generic taxa: {0}".format(", ".join(generics2)),
            "Candidate total: {0} votes".format(self.total_votes),
            "WINNER: {0}".format(self.winning_candidate.to_string()),
            "",
        ]
        message = "\n".join(parts)
        logger = logging.getLogger("brocc.votes")
        logger.debug(message)


class NoAssignment(object):
    is_valid_assignment = False

    """Null object representing no assignment, with message."""
    def __init__(self, query_id, message, rank=None, candidates=None, generics=None):
        self.query_id = query_id
        self.message = message
        self.rank = rank
        self.candidates = candidates
        self.generics = generics

    @property
    def total_votes(self):
        return sum(c.votes for c in self.candidates)

    @property
    def num_generic(self):
        return sum(n for n, taxon in self.generics)

    def format_for_standard_taxonomy(self):
        return "%s\t%s\n" % (self.query_id, self.message)

    format_for_log = format_for_standard_taxonomy

    def log_details(self):
        candidates2 = [c.to_string() for c in self.candidates if c.votes > 1]
        generics2 = ["{0} ({1} votes)".format(k, v) for v, k in self.generics]
        parts = [
            "",
            "** {0} voting for {1} **".format(self.query_id, self.rank),
            "Candidates (>1 vote): {0}".format(", ".join(candidates2)),
            "Generic taxa: {0}".format(", ".join(generics2)),
            "Candidate total: {0} votes".format(self.total_votes),
            self.message,
            "",
        ]
        message = "\n".join(parts)
        logger = logging.getLogger("brocc.votes")
        logger.info(message)

class Assigner(object):
    ranks = [
        "species", "genus", "family", "order",
        "class", "phylum", "kingdom", "superkingdom",
        ]

    def __init__(self, min_cover, species_min_id, genus_min_id, min_id,
                 consensus_thresholds, min_winning_votes, taxa_db):
        self.min_cover = min_cover
        self.rank_min_ids = [
            species_min_id, genus_min_id, min_id, min_id,
            min_id, min_id, min_id, min_id,
            ]
        self.min_id = min_id
        self.consensus_thresholds = consensus_thresholds
        self.min_winning_votes = min_winning_votes
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
            message = "Abundance of low coverage hits: possible chimera"
            return NoAssignment(name, message)
        if not hits_to_keep:
            message = "All BLAST hits were filtered for low quality."
            return NoAssignment(name, message)
        return self.vote(name, seq, hits_to_keep)

    def _retrieve_lineage(self, hit):
        taxid = self.taxa_db.get_taxon_id(hit.accession)
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
            a = self.vote_at_rank(name, rank, hits_lineage)
            a.log_details()
            if a.is_valid_assignment:
                return a
        a.message = "Could not find consensus at domain level. No classification."
        return a

    def vote_at_rank(self, query_id, rank, db_hits):
        '''Votes at a given rank of the taxonomy.'''
        rank_idx = self.ranks.index(rank)
        min_pct_id = self.rank_min_ids[rank_idx]
        consensus_threshold = self.consensus_thresholds[rank_idx]

        # We need to make a distinction between three types of
        # assignment candidates as we tally the votes:
        #
        # 1. Normal taxa (e.g. Debaryomycetaceae)
        # 2. Placeholders, created to fill a rank (e.g. Mortierellales (class))
        # 3. Generic taxa (e.g. uncultured Ascomycota)

        candidates = dict()
        num_generic_votes = 0
        generics = collections.defaultdict(int)
        for hit, lineage in db_hits:
            if hit.pct_id <= min_pct_id:
                continue
            taxon = lineage.get_taxon(rank)
            if taxon is None:
                continue
            if lineage.is_generic(rank):
                generics[taxon] += 1
                num_generic_votes += 1
                continue
            if taxon not in candidates:
                candidates[taxon] = AssignmentCandidate(lineage, rank)
            candidates[taxon].votes += 1

        sorted_candidates = list(
            sorted(candidates.values(), reverse=True, key=lambda c: c.votes))
        total_candidate_votes = sum(c.votes for c in sorted_candidates)
        votes_needed_to_win = total_candidate_votes * consensus_threshold

        sorted_generics = list(
            reversed(sorted((v, k) for k, v in generics.items())))

        if len(candidates) == 0:
            return NoAssignment(
                query_id, "No candidates", rank=rank,
                candidates=sorted_candidates, generics=sorted_generics)

        # The generic taxa shouldn't count towards the vote totals.

        leading_candidate = sorted_candidates[0]

        # Placeholder taxa should count towards the total.  However,
        # if a placeholder taxon wins the vote, we'd like to return
        # the normal taxon from which the placeholder is derived.
        # Thus, if a placeholder wins, we return a null result and
        # kick things up to the next rank.
        if leading_candidate.is_placeholder():
            return NoAssignment(
                query_id, "Placeholder taxon", rank=rank,
                candidates=sorted_candidates, generics=sorted_generics)


        if leading_candidate.votes < votes_needed_to_win:
            message = (
                "No consensus: Leading candidate needed {0} votes to win, "
                "had {1}".format(votes_needed_to_win, leading_candidate.votes))
            return NoAssignment(
                query_id, message, rank=rank,
                candidates=sorted_candidates, generics=sorted_generics)

        if leading_candidate.votes < self.min_winning_votes:
            message = (
                "Absolute number of votes too small for leading candidate: "
                "{0} votes observed, need {1} at a minimum".format(
                    leading_candidate.votes, self.min_winning_votes))
            return NoAssignment(
                query_id, message, rank=rank,
                candidates=sorted_candidates, generics=sorted_generics)

        return Assignment(
            query_id, leading_candidate, rank=rank,
            candidates=sorted_candidates, generics=sorted_generics)
