import random

from matplotlib.colors import ListedColormap


def hex_to_cmap(num_to_generate):

    hexes = generate_hex(num_to_generate)
    return ListedColormap(hexes)

def float_to_hex(c, cmap):
    """
    Given a float value, return the hex color representation.

    Parameters
    ----------
    c : float

    Returns
    -------
    hex string

    """
    c = cmap(c)
    return '#%02x%02x%02x' % (c[0] * 255, c[1] * 255, c[2] * 255)

def generate_hex(num_to_generate):
    """

    Parameters
    ----------
    num_to_generate : int
        number of hex to generate

    Returns
    -------
    list of hex colors
    """
    hexes = []
    for i in range(0,num_to_generate):
        hexes.append(
            '#' + ''.join(
                [random.choice('0123456789ABCDEF') for x in range(6)]
            )
        )
    return hexes


def rgb_to_hex(rgb):
    """
    Returns list of html hex codes for each rgb code in rgb list

    Parameters
    ----------
    rgb : list
        list of RGB values

    Returns
    -------
    list of html hex values
    """
    hexcolors = [
        '#%02x%02x%02x' % (c[0] * 255, c[1] * 255, c[2] * 255) for c in rgb
        ]
    return hexcolors


def expr_to_hex(expr, cmap='Purples', is_norm=True):
    """
    From a list of expression values, return a list of html hex values.

    Parameters
    ----------
    expr : list
        list of positive expression values
    cmap : basestring
        colormap string
    is_norm : Boolean
        scales to max of 1 if true

    Returns
    -------
    list of hex corresponding to expression values
    """
    rgbs = expr_to_rgb(expr, cmap, is_norm)
    return rgb_to_hex(rgbs)


def expr_series_to_hex(expr, cmap='Purples', is_norm=True):
    """
    Takes a series of float values and returns a series of corresponding hex
    according to the cmap specified.

    Parameters
    ----------
    expr : pandas.Series
    cmap : basestring
    is_norm : normalize 0 to 1

    Returns
    -------
    pandas.Series

    """
    if is_norm:
        expr = expr + 1
        normed = expr / expr.max()
    else:
        normed = expr
    normed = normed - 0.1
    return normed.apply(
        lambda x: float_to_hex(x, cmap)
    )

def expr_to_rgb(expr, cmap='Purples', is_norm=True):
    """
    From a list of expression values, return a list of rgb values.

    Parameters
    ----------
    expr : list
        list of positive expression values
    cmap : basestring
        colormap string
    is_norm : Boolean
        scales to max of 1 if true

    Returns
    -------
    list of rgb corresponding to expression values
    """
    if is_norm:
        norm = [float(i) / max(expr) for i in expr]  # normalize expr val [0-1]
    else:
        norm = expr

    rgbs = [cmap(color) for color in norm]  # add color
    return rgbs


def color_by_condition(df, col_string):
    """
    Takes df of conditions (rows of samples, cols of conditions)
    and returns a list of colors matching each distinct condition.

    Parameters
    ----------
    df : pandas.DataFrame
        table of samples as rows, conditions as columns
    col_string : basestring
        column for which to map colors to
    Returns
    -------
    dictionary of {sample:{color:COLOR, condition:CONDITION}, } for each sample
    for each condition in df[col_string]
    """
    try:
        max_conditions = set(df[col_string])
    except KeyError:
        print("{} does not exist as key in dataframe".format(col_string))
        print("reverting to {}".format(df.columns[0]))
        max_conditions = set(df[df.columns[0]])
    colormap = {}
    c = 0
    for condition in max_conditions:
        colormap[condition] = {
            'color':c,
            'condition':condition
        }
        c+=1
    colormapped = df[col_string].apply(lambda x: colormap[x])
    return dict(colormapped)
