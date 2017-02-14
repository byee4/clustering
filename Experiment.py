import ExpressionTable
import color_helpers as ch
from ExpressionTable import *


class Experiment():
    """
    Contains expression (ExpressionTable) and metadata (DataFrame) info
    for a cluster analysis. The table contains expression data in the form of
    counts, featurecounts, or expression for each attribute, for each sample.
    The metadata info contains categorical information about each sample. For
    example, ExpressionTable(iris.txt) contains expression for each attribute,
    while the metadata(iris.names) describes each sample name.
    """
    def __init__(self, counts_file, conditions_file=None, conditions_col=None,
                 gene_id=None, is_featurecounts=False):

        if is_featurecounts:
            self.counts = FeatureCountsTable(counts_file)
        else:
            self.counts = ExpressionTable(counts_file)

        self.source = conditions_file
        self.gene_of_interest = gene_id
        self.condition_of_interest = conditions_col

        self.metadata = self.set_metadata(
            conditions_file,
            conditions_col,
            gene_id,
        )


    def set_metadata(self, conditions_file, conditions_col, gene_id):
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
        gene_id : basestring

        Returns
        -------

        """
        if conditions_file is not None and conditions_col is not None:
            return self._generate_metadata_from_conditions()
        elif gene_id is not None:
            return self._generate_metadata_from_gene_expression()
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
        colors = {}
        for col in self.counts.data.columns:
            colors[col] = {'color': 'blue', 'condition': 'condition'}
        return pd.DataFrame(colors).T

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
        colors = ch.color_by_condition(conditions_df, self.condition_of_interest)
        # self.cmap = ch.hex_to_cmap(conditions_df.shape[0]) # this can be done better.
        return pd.DataFrame(colors).T

    def _generate_metadata_from_gene_expression(self):
        """
        given a gene_id of interest, adds colors based on the (LOG2) level of
        expression for that gene, if that gene is described in the counts file.
        Otherwise, returns generic metadata

        @see _generate_metadata_from_nothing()

        Returns
        -------
        pandas.DataFrame containing samples, color, condition information.

        """
        colors = {}
        if self.gene_of_interest in self.counts.data.index:
            expr = self.counts.data.loc[self.gene_of_interest]
            expr = np.log2(expr+1)
            for key, value in expr.iteritems():
                colors[key] = {'color': value,
                               'condition': 'expression'}
            return pd.DataFrame(colors).T

        else:
            print("warning, gene not found in table.")
            return self._generate_metadata_from_nothing()

    def recolor(self, gene_id):
        self.gene_of_interest = gene_id
        self.metadata = self._generate_metadata_from_gene_expression()
