import numpy as np
import pandas as pd
from sklearn import datasets
from sklearn.cluster import KMeans

class ExpressionTable:

    def __init__(self, data_file):
        """

        Parameters
        ----------
        data_file : basestring
            data file containing index in the first column, expression data
            in subsequent columns with header information in the first row.
        lengths_file : basestring

        """
        data = pd.read_table(data_file, index_col=0)

        self._test_iris()
        self.data = data
        self._pseudocount = 0
        self._is_log2 = False
        self._num_samples = self.data.shape[0]

    def _test_iris(self):
        iris = datasets.load_iris()
        petal_data = iris.data[:, 2:]  # get only petal features, which are the third and fourth values in each sample

        # perform k-means analysis on iris data

        # there are only 3 iris flower groups: 'setosa', 'versicolor', 'virginica'

        kmean = KMeans(n_clusters=3)  # n_clusters asks for only 3 groupings
        kmean.fit(petal_data)
        # print(kmean.labels_)
        # pd.DataFrame(kmean.fit(petal_data)).to_csv('/home/bay001/projects/codebase/clustering/examples/iris_example.test.txt')

    def as_log2(self, pseudocount=1):
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
