import pandas as pd

import ExpressionTable
import helper_functions as h
from ExpressionTable import ExpressionTable


class Experiment:
    """
    Contains expression (ExpressionTable) and metadata (DataFrame) info
    for a cluster analysis. The table contains expression data in the form of
    counts, featurecounts, or expression for each attribute, for each sample.
    The metadata info contains categorical information about each sample. For
    example, ExpressionTable(iris.txt) contains expression for each attribute,
    while the metadata(iris.names) describes each sample name.
    """
    def __init__(self, matrix_file, conditions_file=None, conditions_col=None):

        self.counts = ExpressionTable(matrix_file)

        self.source = conditions_file
        self.condition_of_interest = conditions_col

        self.metadata = self.set_metadata(
            conditions_file,
            conditions_col,
        )

    def set_metadata(self, conditions_file, conditions_col):
        """
        Generates metadata for gene expression. First checks for a
        conditions_file + conditions_col, then checks for gene_id,
        lastly creates a mock metadata file if no metadata is present.

        Parameters
        ----------
        conditions_file : basestring
            tab separated file containing sample in rows, conditions in cols
        conditions_col : basestring
            a condition

        Returns
        -------

        """
        if conditions_file is not None and conditions_col is not None:
            return self._generate_metadata_from_conditions()
        else:
            return self._generate_metadata_from_nothing()

    def _generate_metadata_from_nothing(self):
        """
        If no conditions_file or condition is set, or if no gene_id is set,
        return a generic metadata file (all points are blue, no conditions
        are specified).

        Returns
        -------
        pandas.DataFrame containing samples, color, condition information.

        """
        marker_map = {}
        for col in self.counts.data.columns:
            marker_map[col] = {'marker': 'o', 'condition': 'condition'}
        return pd.DataFrame(marker_map, index=self.counts.data.index).T

    def _generate_metadata_from_conditions(self):
        """
        given a conditions file, returns metadata dataframe containing
        condition and color information based on the column information
        (one color for every unique condition in self.condition_of_interest).

        Returns
        -------
        pandas.DataFrame containing samples, color, condition information.

        """
        conditions_df = pd.read_table(
            self.source,
            index_col=0
        )
        shapes = h.shape_by_condition(conditions_df, self.condition_of_interest)
        shapes = pd.DataFrame(shapes)

        shapes.columns = [str(c) for c in shapes.columns]

        shapes = shapes.T
        # self.cmap = ch.hex_to_cmap(conditions_df.shape[0]) # this can be done better.
        return shapes

    def update_metadata(self, df):
        self.metadata = pd.merge(  # update metadata
            self.metadata,
            df,
            how='left',
            left_index=True,
            right_index=True
        )

    def recolor(self, gene_id):
        self.gene_of_interest = gene_id
        self.metadata = self._generate_metadata_from_gene_expression()
