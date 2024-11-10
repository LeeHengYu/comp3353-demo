class MotifFinding:
    """Class for deterministic motif finding."""

    def __init__(self, size=8, seqs=None):
        self.motif_size = size
        self.seqs = seqs if seqs is not None else []
        if self.seqs:
            self.alphabet = "ACGT"  # Assuming DNA sequences only

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, n):
        return self.seqs[n]

    def seq_size(self, i):
        return len(self.seqs[i])

    def create_motif_from_indexes(self, indexes):
        res = [[0] * self.motif_size for _ in range(len(self.alphabet))]

        for i, ind in enumerate(indexes):
            subseq = self.seqs[i][ind:ind + self.motif_size]
            for j in range(self.motif_size):
                for k in range(len(self.alphabet)):
                    if subseq[j] == self.alphabet[k]:
                        res[k][j] += 1
        return res

    def score(self, s):
        score = 0
        mat = self.create_motif_from_indexes(s)
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1, len(mat)):
                if mat[i][j] > maxcol:
                    maxcol = mat[i][j]
            score += maxcol
        return score

    def next_solution(self, s):
        next_sol = [0] * len(s)
        pos = len(s) - 1
        while pos >= 0 and s[pos] == self.seq_size(pos) - self.motif_size:
            pos -= 1
        if pos < 0:
            return None
        else:
            for i in range(pos):
                next_sol[i] = s[i]
            next_sol[pos] = s[pos] + 1
            for i in range(pos + 1, len(s)):
                next_sol[i] = 0
        return next_sol

    def exhaustive_search(self):
        best_score = -1
        res = []
        s = [0] * len(self.seqs)
        while s is not None:
            sc = self.score(s)
            if sc > best_score:
                best_score = sc
                res = s.copy()
            s = self.next_solution(s)
        return res, best_score

    def heuristic_consensus(self):
        """
        Heuristic approach to find the consensus motif.
        
        Returns:
            list: Starting indices of the motif in each sequence.
        """
        res = [0] * len(self.seqs)
        max_score = -1
        partial = [0, 0]

        for i in range(self.seq_size(0) - self.motif_size + 1):
            for j in range(self.seq_size(1) - self.motif_size + 1):
                partial[0] = i
                partial[1] = j
                sc = self.score(partial)
                if sc > max_score:
                    max_score = sc
                    res[0] = i
                    res[1] = j

        for k in range(2, len(self.seqs)):
            partial = [0] * (k + 1)
            for j in range(k):
                partial[j] = res[j]
            max_score = -1

            for i in range(self.seq_size(k) - self.motif_size + 1):
                partial[k] = i
                sc = self.score(partial)
                if sc > max_score:
                    max_score = sc
                    res[k] = i

        return res, max_score
    
if __name__ == "__main__":
    sequences = [
        "ATCGATCGA",
        "GCTCGATCG",
        "TATCGTATC",
        "CGATCGTCA",
    ]
    finder = MotifFinding(size=6, seqs=sequences)

    bf, bf_score = finder.exhaustive_search()
    heu, heu_score = finder.heuristic_consensus()

    print(f"Best motif positions: {bf} with score {bf_score}")
    print(f"Best motif positions: {heu} with score {heu_score}")
