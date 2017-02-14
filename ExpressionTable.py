import numpy as np
import pandas as pd


class ExpressionTable():

    def __init__(self, data_file, lengths_file = None):
        """

        Parameters
        ----------
        data_file : basestring
            data file containing index in the first column, expression data
            in subsequent columns with header information in the first row.
        lengths_file : basestring

        """
        data = pd.read_table(data_file, index_col=0)
        self.data = data

        self._pseudocount = 0
        self._is_log2 = False
        self._is_rpkm = False
        self._num_samples = self.data.shape[0]

        if lengths_file is not None:
            self.length = self.set_lengths(lengths_file)
        else:
            self.length = None

    def as_log2(self, pseudocount=0):
        """
        log2 transforms self.data

        Parameters
        ----------
        pseudocount : int
            number of fake counts to add to prevent log2(0) inf problem

        Returns
        -------

        """
        self._is_log2 = True
        self._pseudocount = pseudocount
        self.data = np.log2(self.data + pseudocount)

    def subset(self, subset_file):
        """
        removes any row (gene) that is NOT contained in the subset file

        Parameters
        ----------
        subset_file : basestring
            line-delimited file name of gene ids to subset the raw matrix on.

        Returns
        -------

        """
        attributes = [attr.strip() for attr in open(subset_file, 'r')]
        self.data = self.data.loc[attributes, :].dropna(axis=0)

    def min_row_sum_cutoff(self, min_expr_sum=0):
        """
        removes any row (gene) that does NOT meet the minimum row sum
        of read counts (or RPKM if flagged) requirements.

        Parameters
        ----------
        min_expr_sum : int
            sum across all samples to filter

        Returns
        -------

        """
        self.data = self.data[self.data.sum(axis=1) >= min_expr_sum]
        if self.length is not None:
            self.length = self.length.ix[self.data.index]

    def as_rpkm(self):
        """
        Transforms data into RPKM

        Parameters
        ----------
        df : pandas.DataFrame
            featureCounts table formatted by default as having the following
            columns preceding expression data:

            Geneid	Chr	Start	End	Strand	Length  expression**

        Returns
        -------
        """

        counts = self.data
        lengths = self.length
        mapped_reads = counts.sum()
        self.is_rpkm = True
        self.data = (
            counts * pow(10, 9)
        ).div(
            mapped_reads,
            axis=1
        ).div(
            lengths,
            axis=0
        )

    def set_lengths(self, lengths_file):
        """

        Parameters
        ----------
        lengths_file : basestring
            tab delimited file containing the gene name and corresponding
            length. Must contain header information.
        Returns
        -------

        """
        self.lengths = pd.read_table(lengths_file,index_col=0,names=['Length'])

        # check if all lengths and expr counts are there, otherwise warn
        if set(lengths.index).intersection(set(self.data.index)
                                           ) != self.num_samples:
            print("Warning: Length annotations and gene expr dont match." +
                  "Taking intersection of {} length features and " +
                  "{} genes".format(
                      len(self.lengths),
                      self.num_samples)
                  )
            self.data = self.data[self.lengths.index].dropna()
            self.lengths = self.lengths.ix[
                set(self.lengths.index).intersection(set(self.data.index))
            ][self.lengths.columns[0]]


class FeatureCountsTable(ExpressionTable):
    """
    This class uses featurecounts counts.txt to populate gene and gene len
    info.
    """
    def __init__(self, counts_file):
        """

        Parameters
        ----------
        counts_file : basestring
            featureCounts counts.txt
        """
        counts = pd.read_table(counts_file, index_col=0, skiprows=1)
        self.data = counts.ix[:,5:]
        self.length = counts['Length']

        self._pseudocount = 0
        self._is_rpkm = False
