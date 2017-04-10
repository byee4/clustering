import matplotlib

matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns
from bokeh.models import ColumnDataSource

import helper_functions as h

__all__ = []
__version__ = 0.1
__date__ = '2017-2-14'
__updated__ = '2017-2-14'


class _Plotter():

    def __init__(self, expt, n_clusters, seed, cmap = 'Purples'):
        """

        Parameters
        ----------
        data : pandas.DataFrame
            A table of gene expression in the format (genes, samples)
        colors : pandas.Dataframe
            A table defining the color and condition for each sample name

        """
        self.expt = expt
        self.seed = seed
        self.n_clusters = n_clusters
        self.cmap = plt.get_cmap(cmap)
        self.data = expt.counts.data
        self.kmeans, self.kclust = self._fit_transform()

    def _merge_metadata_and_cluster_data(self):
        """
        Merges the metadata (labels and such) with the cluster data

        Returns
        -------
        merged : pandas.DataFrame

        """
        merged = pd.merge(
            self.expt.metadata,
            self.kclust,
            how='left',
            left_index=True,
            right_index=True
        )
        return merged

    def _fit_transform(self):
        """
        Transforms the data into k clusters and updates expt metadata.

        Returns
        -------
        transf : pandas.DataFrame
            table containing kmeans values
        """
        if self.seed != 0:
            k = KMeans(n_clusters=self.n_clusters, random_state=self.seed)
        else:
            k = KMeans(n_clusters=self.n_clusters)
        transf = k.fit_transform(self.data.T)
        # fit = k.fit_predict(self.data.T)
        # print(fit[:100]) # TODO change shape based on classification.
        transf = pd.DataFrame(transf, index=self.data.columns)

        labels = pd.DataFrame(
            k.labels_, columns=['label'], index=self.data.columns
        )

        self.expt.update_metadata(labels)
        return k, transf

    def _matplotlib(self, ax=None):
        """

        Parameters
        ----------
        ax : matplotlib.axes._subplots.AxesSubplot
            subplot axes

        Returns
        -------

        """
        if ax is None:
            ax = plt.gca()

        df = self._merge_metadata_and_cluster_data()

        colors = sns.color_palette("hls", len(set(df['label'])))

        for _, row in df.iterrows():
            ax.scatter(
                row[0], row[1],
                label=row['condition'],
                color=colors[row['label']],
                marker=row['marker']
            )

    def plot(self, bokeh=False, ax=None):
        """

        Parameters
        ----------
        bokeh : Boolean
            True if plotting bokeh figure, else matplotlib axes
        ax : matplotlib.axes._subplots.AxesSubplot or bokeh.plotting.figure.Figure

        Returns
        -------

        """
        if bokeh:
            self._bokeh(ax)
        else:
            self._matplotlib(ax)

    def update_cmap(self):
        pass


def kmeansplot(expt, n_clusters, cmap, seed, ax=None, bokeh=False):
    """

    Parameters
    ----------
    expt : Experiment
        Object defining the expression data and conditions for samples.
    cmap : basestring
        colormap string
    ax : matplotlib.axes._subplots.AxesSubplot or bokeh.plotting.figure.Figure
    bokeh : Boolean
        True if plotting bokeh figure, else matplotlib axes

    Returns
    -------
    _PCAPlotter object

    """

    plotter = _Plotter(expt, n_clusters, seed, cmap)
    plotter.plot(bokeh=bokeh, ax=ax)
    return plotter