from __future__ import division
import collections

from brocclib.taxonomy import Lineage, NoLineage, is_generic

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
    def all_taxa(self):
        return self.lineage.get_all_taxa(self.rank)

    @property
    def standard_taxa(self):
        return self.lineage.get_standard_taxa(self.rank)

    def to_string(self):
        return "<AssignmentCandidate {0} {1} {2}>".format(self.rank, self.lineage.get_taxon(self.rank), self.votes)


class Assignment(object):
    is_valid_assignment = True

    def __init__(self, query_id, winning_candidate, total_votes, num_generic):
        self.query_id = query_id
        self.winning_candidate = winning_candidate
        self.total_votes = total_votes
        self.num_generic = num_generic

    def format_for_full_taxonomy(self):
        lineage = ';'.join(self.winning_candidate.all_taxa)
        return "%s\t%s\n" % (self.query_id, lineage)

    def format_for_standard_taxonomy(self):
        lineage = ';'.join(self.winning_candidate.standard_taxa)
        return "%s\t%s\n" % (self.query_id, lineage)

    def format_for_log(self):
        lineage = ';'.join(self.winning_candidate.all_taxa)
        return "%s\t%s\t%s\t%s\t%s\t%s\n" % (
            self.query_id, self.winning_candidate.votes, self.total_votes,
            self.num_generic, self.winning_candidate.rank, lineage)


class NoAssignment(object):
    is_valid_assignment = False

    """Null object representing no assignment, with message."""
    def __init__(self, query_id, message):
        self.query_id = query_id
        self.message = message

    def format_for_full_taxonomy(self):
        return "%s\t%s\n" % (self.query_id, self.message)

    format_for_standard_taxonomy = format_for_full_taxonomy
    format_for_log = format_for_full_taxonomy


class Assigner(object):
    ranks = [
        "species", "genus", "family", "order",
        "class", "phylum", "kingdom", "domain",
        ]

    def __init__(self, min_cover, species_min_id, genus_min_id, min_id,
                 consensus_thresholds, min_votes, taxa_db):
        self.min_cover = min_cover
        self.rank_min_ids = [
            species_min_id, genus_min_id, min_id, min_id,
            min_id, min_id, min_id, min_id,
            ]
        self.min_id = min_id
        self.consensus_thresholds = consensus_thresholds
        self.min_votes = min_votes
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
            if a.is_valid_assignment:
                return a
        return NoAssignment(
            name, "Could not find consensus at domain level. No classification.")

    def vote_at_rank(self, query_id, rank, db_hits):
        '''Votes at a given rank of the taxonomy.'''
        # print
        # print "====", query_id, "===="
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
            if is_generic(taxon):
                generics[taxon] += 1
                num_generic_votes += 1
                continue
            if taxon not in candidates:
                candidates[taxon] = AssignmentCandidate(lineage, rank)
            candidates[taxon].votes += 1

        if len(candidates) == 0:
            return NoAssignment(query_id, "No candidates")

        total_candidate_votes = sum(c.votes for c in candidates.values())
        votes_needed_to_win = total_candidate_votes * consensus_threshold

        sorted_candidates = candidates.values()
        sorted_candidates.sort(reverse=True, key=lambda c: c.votes)

        # print [c.to_string() for c in sorted_candidates if c.votes > 0]
        # print "Generics", list(reversed(sorted((v, k) for k, v in generics.iteritems())))
        # print "Votes:", total_candidate_votes, "candidate", num_generic_votes, "generic"
        # print "Need", votes_needed_to_win, "to win", "(threshold", consensus_threshold, ")"

        # The generic taxa shouldn't count towards the vote totals.
        if total_candidate_votes < self.min_votes:
            message = (
                "Not enough votes observed after removing generic taxa: " \
                "{0} candidate votes, {1} generic taxa".format(
                    total_candidate_votes, num_generic_votes))
            return NoAssignment(query_id, message)

        leading_candidate = sorted_candidates.pop(0)

        # Placeholder taxa should count towards the total.  However,
        # if a placeholder taxon wins the vote, we'd like to return
        # the normal taxon from which the placeholder is derived.
        # Thus, if a placeholder wins, we return None and kick things
        # up to the next rank.
        if leading_candidate.is_placeholder():
            return NoAssignment(query_id, "Placeholder taxon")

        if leading_candidate.votes >= votes_needed_to_win:
            return Assignment(query_id, leading_candidate, total_candidate_votes, num_generic_votes)
        else:
            message = "Could not find consensus at {0} level. " \
                  "No classification.".format(rank)
            return NoAssignment(query_id, message)
