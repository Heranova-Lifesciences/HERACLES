from Bio import Align

# Local Smith-Waterman aligner replacing the deprecated Bio.pairwise2 module.
# Parameterisation is identical to the previous pairwise2.align.localms
# (match 2, mismatch -1, gap open -5, gap extend -0.1); scores were verified
# to be numerically identical between the two implementations.
_ALIGNER = Align.PairwiseAligner(
    mode="local",
    match_score=2,
    mismatch_score=-1,
    open_gap_score=-5,
    extend_gap_score=-0.1,
)


def _render_with_ends(coreA, coreB, sA, eA, sB, eB, seqA, seqB):
    """Embed the local-alignment core between the unaligned ends.

    pairwise2 rendered local alignments over the full length of both
    sequences: every unaligned residue occupies a column opposite a gap
    in the other string.  The predicates depend on that rendering (the
    cleaned strings are always the original sequences), so it is kept.
    """
    preA, preB = seqA[:sA], seqB[:sB]
    sufA, sufB = seqA[eA:], seqB[eB:]
    lead = max(len(preA), len(preB))
    tail = max(len(sufA), len(sufB))
    leadA = preA.rjust(lead, "-")
    leadB = preB.rjust(lead, "-")
    tailA = sufA.ljust(tail, "-")
    tailB = sufB.ljust(tail, "-")
    return leadA + coreA + tailA, leadB + coreB + tailB


def _local_alignments(seqA, seqB):
    """Return all optimal local alignments as (gappedA, gappedB, score).

    The gapped strings reproduce the rendering that Bio.pairwise2 used for
    local alignments (see _render_with_ends).  Scores are the true local
    (Smith-Waterman) optima, numerically identical to pairwise2's.
    """
    out = []
    for aln in _ALIGNER.align(seqA, seqB):
        blocks = aln.aligned
        if blocks[0].size == 0:
            out.append(("-" * len(seqA), "-" * len(seqB), aln.score))
            continue
        sA, eA = int(blocks[0][0][0]), int(blocks[0][-1][1])
        sB, eB = int(blocks[1][0][0]), int(blocks[1][-1][1])
        ga, gb = _render_with_ends(aln[0], aln[1], sA, eA, sB, eB,
                                   seqA, seqB)
        out.append((ga, gb, aln.score))
    return out

def has_only_terminal_gaps(seqA, seqB):
    """
    Returns True if gaps occur only at the ends (sticky ends),
    and there are no internal gaps in the aligned region.
    """
    # find aligned region (where both have bases)
    aligned_positions = [
        i for i, (a, b) in enumerate(zip(seqA, seqB))
        if a != "-" and b != "-"
    ]

    if not aligned_positions:
        return False

    start = aligned_positions[0]
    end = aligned_positions[-1]

    # check no gaps inside aligned region
    for i in range(start, end + 1):
        if seqA[i] == "-" or seqB[i] == "-":
            return False

    return True

def valid_terminal_alignment(seqA, seqB):
    alignments = _local_alignments(seqA, seqB)

    for ga, gb, _ in alignments:
        if has_only_terminal_gaps(ga, gb):
            return True

    return False


def valid_terminal_alignment2(seqA, seqB, min_similarity=0.9):
    """
    Return True if seqA and seqB have a valid terminal alignment
    and the percentage similarity is >= min_similarity.
    
    min_similarity: fraction between 0 and 1
    """
    alignments = _local_alignments(seqA, seqB)

    for ga, gb, _ in alignments:
        if has_only_terminal_gaps(ga, gb):
            # Count aligned positions ignoring gaps
            matches = 0
            aligned_positions = len(ga)
            for x, y in zip(ga, gb):
                if x != '-' and y != '-':
                    if x == y:
                        matches += 1

            if aligned_positions == 0:
                continue

            similarity = matches / aligned_positions  # fraction
            if similarity >= min_similarity:
                return True  # alignment passes

    return False

def has_min_overlap(seqA, seqB, min_overlap=17):
    """
    Return True if seqA and seqB share at least min_overlap
    contiguous nucleotides anywhere.
    """
    lenA, lenB = len(seqA), len(seqB)
    # Try aligning seqA against every possible substring of seqB of length >= min_overlap
    for startB in range(lenB - min_overlap + 1):
        for startA in range(lenA - min_overlap + 1):
            # Compare min_overlap length
            windowA = seqA[startA:startA+min_overlap]
            windowB = seqB[startB:startB+min_overlap]
            if windowA == windowB:
                return True
    return False

def valid_terminal_alignment_tRF3(seqA, seqB, min_similarity=0.94):
    """
    Return True if seqA and seqB have a valid terminal alignment
    and the percentage similarity is >= min_similarity.
    
    min_similarity: fraction between 0 and 1
    """
    alignments = _local_alignments(seqA, seqB)

    for ga, gb, _ in alignments:
        if has_only_terminal_gaps(ga, gb):
            # Count aligned positions ignoring gaps
            matches = 0
            total_aligned_positions = max(len(ga), len(gb))
            for x, y in zip(ga, gb):
                if x != '-' and y != '-':
                    if x == y:
                        matches += 1

            if total_aligned_positions == 0:
                continue

            similarity = matches / total_aligned_positions  # fraction

            seqA_clean = ga.replace("-", "")
            seqB_clean = gb.replace("-", "")
            nt_aligned_positions = min(len(seqA_clean), len(seqB_clean))

            if similarity >= min_similarity and (seqA_clean[-nt_aligned_positions:] == seqB_clean[-nt_aligned_positions:]):
                return True

    return False

def valid_terminal_alignment_tRF5(seqA, seqB, min_similarity=0.94):
    """
    Return True if seqA and seqB have a valid terminal alignment
    and the percentage similarity is >= min_similarity.
    
    min_similarity: fraction between 0 and 1
    """
    
    alignments = _local_alignments(seqA, seqB)
    
    for ga, gb, _ in alignments:
        if has_only_terminal_gaps(ga, gb):
            # Count aligned positions ignoring gaps
            matches = 0
            total_aligned_positions = max(len(ga), len(gb))
            for x, y in zip(ga, gb):
                if x != '-' and y != '-':
                    if x == y:
                        matches += 1

            if total_aligned_positions == 0:
                continue

            similarity = matches / total_aligned_positions  # fraction

            seqA_clean = ga.replace("-", "")
            seqB_clean = gb.replace("-", "")
            nt_aligned_positions = min(len(seqA_clean), len(seqB_clean))

            if similarity >= min_similarity and (seqA_clean[0:nt_aligned_positions] == seqB_clean[0:nt_aligned_positions]):
                return True

    return False

def kmers(seq, k=10):
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}

def candidate_sequences_tRF5(seq, tRF5_index, seq_dict_tRF5):
    prefix = seq[:7]
    return {other for other in tRF5_index[prefix] if other in seq_dict_tRF5 and other != seq}

def collapse_cluster_tRF5(nodes, tRF5_G):
    rep = max(nodes, key=lambda x: tRF5_G.nodes[x]['count'])
    total_count = sum(tRF5_G.nodes[n]['count'] for n in nodes)
    return rep, total_count

def candidate_sequences_tRF3(seq, tRF3_index, seq_dict_tRF3):
    suffix = seq[-7:]
    return {other for other in tRF3_index[suffix] if other in seq_dict_tRF3 and other != seq}

def collapse_cluster_tRF3(nodes, tRF3_G):
    rep = max(nodes, key=lambda x: tRF3_G.nodes[x]['count'])
    total_count = sum(tRF3_G.nodes[n]['count'] for n in nodes)
    return rep, total_count

def valid_terminal_alignment_tRFi(seqA, seqB, min_similarity=0.94):
    """
    Return True if seqA and seqB have a valid terminal alignment
    and the percentage similarity is >= min_similarity.
    
    min_similarity: fraction between 0 and 1
    """
    alignments = _local_alignments(seqA, seqB)

    for ga, gb, _ in alignments:
        if has_only_terminal_gaps(ga, gb):
            # Count aligned positions ignoring gaps
            matches = 0
            total_aligned_positions = max(len(ga), len(gb))
            for x, y in zip(ga, gb):
                if x != '-' and y != '-':
                    if x == y:
                        matches += 1

            if total_aligned_positions == 0:
                continue

            similarity = matches / total_aligned_positions  # fraction

            start = None
            end = None

            for i, (c, d) in enumerate(zip(ga, gb)):
                if c != '-' and d != '-':
                    if start is None:
                        start = i
                    end = i

            if start is None:
                continue
            
            if similarity >= min_similarity and (ga[start:end+1] == gb[start:end+1]):
                return True

    return False

def collapse_cluster_tRFi(nodes, G):
    rep = max(nodes, key=lambda x: G.nodes[x]['count'])
    total_count = sum(G.nodes[n]['count'] for n in nodes)
    return rep, total_count


def adaptive_subset(seq, min_k=7, max_k=12):
    seq_len = len(seq)
    if seq_len < min_k:
        return seq
    k = min(max_k, max(min_k, seq_len // 3))
    start = (seq_len - k) // 2
    return seq[start:start + k]


def align_to_ref(query, matches_dict, original_tsRNA):
    results = []
    for target_seq, target_count in matches_dict.items():
        alignments = _local_alignments(query, target_seq)
        for ga, gb, score in alignments:
            matches_count = sum(1 for a, b in zip(ga, gb) if a == b and a != '-')
            aligned_len = max(len(ga), len(gb))
            similarity = matches_count / aligned_len if aligned_len > 0 else 0

            results.append({
                'original_tsRNA': original_tsRNA,
                'query_seq': query,
                'aligned_tsRNA': target_seq,
                'aligned_tsRNA_similar_cleaned': target_seq,
                'target_count': target_count,
                'similarity': round(similarity, 4),
                'alignment_score': int(score)
            })
            break

    results.sort(key=lambda x: x['alignment_score'], reverse=True)
    return results