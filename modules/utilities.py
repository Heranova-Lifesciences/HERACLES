from Bio import pairwise2

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
    alignments = pairwise2.align.localms(
        seqA, seqB,
        2,    # match
        -1,   # mismatch
        -5,   # gap open
        -0.1  # gap extend
    )

    for a in alignments:
        if has_only_terminal_gaps(a.seqA, a.seqB):
            return True

    return False


def valid_terminal_alignment2(seqA, seqB, min_similarity=0.9):
    """
    Return True if seqA and seqB have a valid terminal alignment
    and the percentage similarity is >= min_similarity.
    
    min_similarity: fraction between 0 and 1
    """
    alignments = pairwise2.align.localms(
        seqA, seqB,
        2,    # match
        -1,   # mismatch
        -5,   # gap open
        -0.1  # gap extend
    )

    for a in alignments:
        if has_only_terminal_gaps(a.seqA, a.seqB):
            # Count aligned positions ignoring gaps
            matches = 0
            aligned_positions = len(a.seqA)
            for x, y in zip(a.seqA, a.seqB):
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
    alignments = pairwise2.align.localms(
        seqA, seqB,
        2,    # match
        -1,   # mismatch
        -5,   # gap open
        -0.1  # gap extend
    )

    for a in alignments:
        if has_only_terminal_gaps(a.seqA, a.seqB):
            # Count aligned positions ignoring gaps
            matches = 0
            total_aligned_positions = max(len(a.seqA), len(a.seqB))
            for x, y in zip(a.seqA, a.seqB):
                if x != '-' and y != '-':
                    if x == y:
                        matches += 1

            if total_aligned_positions == 0:
                continue

            similarity = matches / total_aligned_positions  # fraction

            nt_aligned_positions = min(len(a.seqA.replace("-", "")), len(a.seqB.replace("-", "")))           
            
            if similarity >= min_similarity and (a.seqA[-nt_aligned_positions:] == a.seqB[-nt_aligned_positions:]):
                return True  # alignment passes
            else:
                return False

def valid_terminal_alignment_tRF5(seqA, seqB, min_similarity=0.94):
    """
    Return True if seqA and seqB have a valid terminal alignment
    and the percentage similarity is >= min_similarity.
    
    min_similarity: fraction between 0 and 1
    """
    
    alignments = pairwise2.align.localms(
        seqA, seqB,
        2,    # match
        -1,   # mismatch
        -5,   # gap open
        -0.1  # gap extend
    )
    
    for a in alignments:
        if has_only_terminal_gaps(a.seqA, a.seqB):
            # Count aligned positions ignoring gaps
            matches = 0
            total_aligned_positions = max(len(a.seqA), len(a.seqB))
            for x, y in zip(a.seqA, a.seqB):
                if x != '-' and y != '-':
                    if x == y:
                        matches += 1

            if total_aligned_positions == 0:
                continue

            similarity = matches / total_aligned_positions  # fraction

            nt_aligned_positions = min(len(a.seqA.replace("-", "")), len(a.seqB.replace("-", "")))           
            
            if similarity >= min_similarity and (a.seqA[0:nt_aligned_positions] == a.seqB[0:nt_aligned_positions]):
                return True  # alignment passes
            else:
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
    alignments = pairwise2.align.localms(
        seqA, seqB,
        2,    # match
        -1,   # mismatch
        -5,   # gap open
        -0.1  # gap extend
    )

    for a in alignments:
        if has_only_terminal_gaps(a.seqA, a.seqB):
            # Count aligned positions ignoring gaps
            matches = 0
            total_aligned_positions = max(len(a.seqA), len(a.seqB))
            for x, y in zip(a.seqA, a.seqB):
                if x != '-' and y != '-':
                    if x == y:
                        matches += 1

            if total_aligned_positions == 0:
                continue

            similarity = matches / total_aligned_positions  # fraction

 #           nt_aligned_positions = min(len(a.seqA.replace("-", "")), len(a.seqB.replace("-", "")))      this should be removed,      
            start = None
            end = None

            for i, (c, d) in enumerate(zip(a.seqA, a.seqB)):
                if c != '-' and d != '-':
                    if start is None:
                        start = i
                    end = i

            if start is None:
                continue
            
            if similarity >= min_similarity and (a.seqA[start:end+1] == a.seqB[start:end+1]):
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
        alignments = pairwise2.align.localms(
            query, target_seq,
            2, -1, -5, -0.1
        )
        for aln in alignments:
            matches_count = sum(1 for a, b in zip(aln.seqA, aln.seqB) if a == b and a != '-')
            aligned_len = max(len(aln.seqA), len(aln.seqB))
            similarity = matches_count / aligned_len if aligned_len > 0 else 0

            results.append({
                'original_tsRNA': original_tsRNA,
                'query_seq': query,
                'aligned_tsRNA': target_seq,
                'aligned_tsRNA_similar_cleaned': target_seq,
                'target_count': target_count,
                'similarity': round(similarity, 4),
                'alignment_score': int(aln.score)
            })
            break

    results.sort(key=lambda x: x['alignment_score'], reverse=True)
    return results