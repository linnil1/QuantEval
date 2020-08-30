import numpy as np


class Match:
    """
    Record the connection between two sequences.

    The connection is based on blast result.

    Attributes
    ----------
    name: str
        The id of this two sequences (Included orientation)
    q_name: str
        The id of this query sequence
    r_name: str
        The id of this reference sequence

    q_depth: array of float
        The alignment identity of each base on query sequence
    r_depth: array of float
        The alignment identity of each base on reference sequence

    q_global_identity: float
        The identity value of query sequences from all fragments between this two sequences
    p_global_identity: float
        The identity value of reference sequences from all fragments between this two sequences

    """
    def __init__(self, row, q_seq, r_seq):
        """
        Parameters
        ----------
        row: pandas.Series
            The table of blast result
        q_seq: pandas.Series
            The row of sequence in Sequences table
        r_seq: pandas.Series
            The row of sequence in Sequences table
        """
        self.q_name = q_seq['index']
        self.r_name = r_seq['index']
        self.q_depth = np.zeros(q_seq['length']).astype(float)
        self.r_depth = np.zeros(r_seq['length']).astype(float)
        self.q_global_identity = 0
        self.r_global_identity = 0
        self.name = Match.getName(row)

    def extend(self, row):
        """
        Recalculate the identity score for each base
        from the row of blast result.

        'q_depth' and 'r_depth' will be updated

        Parameters
        ----------
        row: pandas.dataframe
            The table of blast result
        """
        # Read from blast record
        q_start  = row.q_start
        q_end    = row.q_end
        r_start  = row.r_start
        r_end    = row.r_end
        identity = row.identity

        # make sure python slice will not fail
        if r_start > r_end:
            r_start, r_end = r_end, r_start

        # Update the alignment map
        self.q_depth[q_start-1:q_end] = np.maximum(self.q_depth[q_start-1:q_end], identity)
        self.r_depth[r_start-1:r_end] = np.maximum(self.r_depth[r_start-1:r_end], identity)

    def calculateIdentity(self):
        """
        Recalculate the identity score between two sequences

        'q_global_identity' and 'r_global_identity' will be updated

        Parameters
        ----------
        row: pandas.dataframe
            The table of blast result
        """
        self.q_global_identity = np.mean(self.q_depth)
        self.r_global_identity = np.mean(self.r_depth)

    @classmethod
    def getName(cls, row):
        """
        Get the name of the pair sequences from a record of blast result

        The name format will be {q_name}.{r_name}.{orientation}

        Parameters
        ----------
        row: pandas.Series
            The table of blast result

        Returns
        -------
        name: str
        """
        if row.r_start > row.r_end:
            orientation = 'R'
        else:
            orientation = 'F'
        return f"{row.q_name}.{row.r_name}.{orientation}"

    def __str__(self):
        return f"<Match {self.name}>"

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, target):
        return isinstance(target, Match) and self.name == target.name
