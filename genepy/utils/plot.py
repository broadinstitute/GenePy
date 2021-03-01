# Jeremie Kalfon
# for BroadInsitute
# in 2019

from __future__ import print_function
import venn as pyvenn
from matplotlib import pyplot as plt
from matplotlib import cm
from bokeh.palettes import *
from bokeh.layouts import layout, widgetbox, column, row
from bokeh.models.annotations import LabelSet
from bokeh.models.widgets import TextInput
from bokeh.models import HoverTool, CustomJS, BasicTicker, ColorBar, PrintfTickFormatter
from bokeh.models import ColumnDataSource, LinearColorMapper, LogColorMapper
from bokeh.util.hex import hexbin
from bokeh.transform import linear_cmap
from bokeh.io import show
from bokeh.plotting import *
import bokeh
import colorcet as cc
from PIL import Image, ImageDraw, ImageFont
import seaborn as sns

import pandas as pd
from math import pi
import numpy as np
import pdb

from taigapy import TaigaClient
tc = TaigaClient()


def scatter(data, labels=None, title='scatter plot', showlabels=False, folder='',
            colors=None, xname='', yname="", importance=None, radi=5, alpha=0.8, **kwargs):
    """
    Makes an interactive scatter plot using Bokeh

    Args:
    -----
      data: array-like with shape [N,2]
      labels: list[str] a list of N names for each points
      title: str the plot title
      showlabels: bool if the labels should be always displayed or not (else just on hover)
      colors: list[int] of N integers from 0 up to 256 for the dot's colors
      folder: str of location where to save the plot, won't save if empty
      xname: str the name of the x axes
      yname: str the name of the y axes
      importance: a list[int] of N values to scale the size of the dots and their opacity by
      radi: int the size of the dots
      alpha: float the opacity of the dots
      **kwargs: additional bokeh.figure args

    Returns:
    ------
      the bokeh object
    """
    TOOLS = "hover,crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,save,box_select,lasso_select,"

    col = viridis(len(set(colors))) if colors is not None else [
        '#29788E']  # (viridis1)
    radii = []
    fill_alpha = []
    cols = []
    for i in range(data.shape[0]):
        radii.append(radi if importance is None else radi /
                     2 + importance[i] * 30)
        fill_alpha.append(
            alpha if importance is None else alpha - (0.2 * importance[i]))
        cols.append(col[0] if colors is None else col[int(colors[i])])
    source = ColumnDataSource(data=dict(
        x=data[:, 0],
        y=data[:, 1],
        labels=labels if labels is not None else [''] * len(radii),
        fill_color=cols,
        fill_alpha=fill_alpha,
        radius=radii
    ))
    TOOLTIPS = [
        ("name", "@labels"),
        ("(x,y)", "(@x, @y)"),
    ]
    p = figure(tools=TOOLS, tooltips=TOOLTIPS, title=title)
    p.circle('x', 'y', color='fill_color',
             fill_alpha='fill_alpha',
             line_width=0,
             radius='radius' if radi else None, source=source)
    p.xaxis[0].axis_label = xname
    p.yaxis[0].axis_label = yname
    if showlabels:
        labels = LabelSet(x='x', y='y', text='labels', level='glyph', text_font_size='9pt',
                          x_offset=5, y_offset=5, source=source, render_mode='canvas')
        p.add_layout(labels)
    try:
        show(p)
    except:
        show(p)
    if folder:
      save(p, folder + title.replace(' ', "_") + "_scatter.html")
      #export_png(p, filename=folder + title.replace(' ', "_") + "_scatter.png")
    return p


def bigScatter(data, precomputed=False, logscale=False, features=False,
               title="BigScatter", binsize=0.1, folder="", showpoint=False):
    """
    uses a binning method to display millions of points at the same time and showcase density.

    Does this in an interactive fashion

    Args:
    -----
      data: array like. the array containing the point location x,y or their location and density of bins (q,r,c)
      precomputed: bool whether or not the array has aleady been hexbined
      logscale: bool, whether or not the data is logscaled
      features: list[] if the matrix contains a feature column with feature information (names, values. for each bin/dot)
      title: str the title of the plot
      binsize: float the size of the bins
      showpoint: bool whether or not to display them as points or hexes.
      folder: str of location where to save the plot, won't save if empty

    Returns:
    ------
      the bokeh object
    """
    TOOLS = "wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,save,box_select,lasso_select,"
    names = [("count", "@c")]
    if features:
        names.append(('features', '@features'))
    if precomputed:
        TOOLS = "hover," + TOOLS
    p = figure(title=title, tools=TOOLS, tooltips=names if precomputed else None,
               match_aspect=True, background_fill_color='#440154')
    if precomputed:
        p.hex_tile(q="q", r="r", size=binsize, line_color=None, source=data,
                   hover_color="pink", hover_alpha=0.8,
                   fill_color=linear_cmap('c', 'Viridis256', 0, max(data.c)) if not logscale else {'field': 'c', 'transform': LogColorMapper('Viridis256')})
    else:
        if features:
            print("we cannot yet process features on non precomputed version")
        r, bins = p.hexbin(data[:, 0], data[:, 1], line_color=None, size=binsize,
                           hover_color="pink", hover_alpha=0.8,
                           fill_color=linear_cmap('c', 'Viridis256', 0, None) if not logscale else {'field': 'c', 'transform': LogColorMapper('Viridis256')})
    p.grid.visible = False
    if showpoint:
        p.circle(data[:, 0], data[:, 1], color="white", size=1)

    if not precomputed:
        p.add_tools(HoverTool(
            tooltips=names,
            mode="mouse", point_policy="follow_mouse", renderers=[r] if not precomputed else None
        ))

    try:
        show(p)
    except:
        show(p)
    if folder:
      save(p, folder + title.replace(' ', "_") + "_scatter.html")
      #export_png(p, filename=folder + title.replace(' ', "_") + "_scatter.png")

    return p


def CNV_Map(df, sample_order=[], title="CN heatmaps sorted by SMAD4 loss, pointing VPS4B",
            width=900, height=400, standoff=10, y_label='', marks=[]):
    """
    create an interactive plot suited for visualizing segment level CN data for a set of samples using bokeh

    args:
    ----
      df: df['Sample' 'Start' 'End' 'Segment_Mean' 'size'] the df containing segment level
        copy number (can be subsetted to a specific region or genome-wide)
      sampleorder: list[Sample] <- for all samples present in the df
      title: plot title
      width: int width
      height: int height
      standoff: the space between the plot and the x axis
      y_label: the y axis label
      marks: location of lines at specific loci

    Returns:
    --------
      The bokeh object
    """
    colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2",
              "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
    colors = RdBu[8]
    mapper = LinearColorMapper(
        palette=colors, low=df.Segment_Mean.min(), high=df.Segment_Mean.max())
    if len(sample_order) == 0:
        sample_order = list(set(df.Sample.tolist()))
    TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"
    p = figure(title=title,
               y_range=(df.End.max(), df.Start.min()),
               x_range=sample_order,
               x_axis_location="above", plot_width=width, plot_height=height,
               tools=TOOLS, toolbar_location='below',
               tooltips=[('pos', '@Start, @End'), ('relative CN', '@Sample')])

    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "5pt"
    p.axis.major_label_standoff = standoff
    p.xaxis.major_label_orientation = pi / 3
    pos = 0
    # for i,val in enumerate(historder):
    #    p.rect(x=pos,y=-7,width=len(orderbyhist[val]), height=10, fill_color=small_palettes['Viridis'][4][i])
    #    p.text(x=pos+len(orderbyhist[val])/2, y=-9, text=str(val), text_color="#96deb3")
    #    pos += len(orderbyhist[val])

    p.rect(x="Sample", y="Start", width=0.9, height="size",
           source=df.reset_index(drop=True),
           fill_color={'field': 'Segment_Mean', 'transform': mapper},
           line_color=None)

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                         ticker=BasicTicker(desired_num_ticks=len(colors)),
                         formatter=PrintfTickFormatter(format="%.2f"),
                         label_standoff=6, border_line_color=None, location=(0, 0))
    p.add_layout(color_bar, 'right')
    p.yaxis.axis_label = y_label
    # p.yaxis.major_label_overrides={20:'Centromer'}
    for val in marks:
        hline = Span(location=val, dimension='width',
                     line_color='green', line_width=0.2)
        p.renderers.extend([hline])
    if folder:
      save(p, folder + title.replace(' ', "_") + "_cn_plot.html")
      #export_png(p, filename=folder + title.replace(' ', "_") + "_cn_plot.png")
    show(p)      # show the plot
    return p


def volcano(data, folder='', tohighlight=None, tooltips=[('gene', '@gene_id')],
            title="volcano plot", xlabel='log-fold change', ylabel='-log(Q)', maxvalue=100,
            searchbox=False, logfoldtohighlight=0.15, pvaltohighlight=0.1, showlabels=False):
    """
    Make an interactive volcano plot from Differential Expression analysis tools outputs

    Args:
    -----
      data: a df with rows genes and cols [log2FoldChange, pvalue, gene_id]
      folder: str of location where to save the plot, won't save if empty
      tohighlight: list[str] of genes to highlight in the plot
      tooltips: list[tuples(str,str)] if user wants tot specify another bokeh tooltip
      title: str plot title
      xlabel: str if user wants to specify the title of the x axis
      ylabel: str if user wants tot specify the title of the y axis
      maxvalue: float the max -log2(pvalue authorized usefull when managing inf vals)
      searchbox: bool whether or not to add a searchBox to interactively highlight genes
      logfoldtohighlight: float min logfoldchange when to diplay points
      pvaltohighlight: float min pvalue when to diplay points
      showlabels: bool whether or not to show a text above each datapoint with its label information

    Returns:
    --------
      The bokeh object
    """
    # pdb.set_trace()
    to_plot_not, to_plot_yes = selector(data, tohighlight if tohighlight is not None else [
    ], logfoldtohighlight, pvaltohighlight)
    hover = bokeh.models.HoverTool(tooltips=tooltips,
                                   names=['circles'])

    # Create figure
    p = bokeh.plotting.figure(title=title, plot_width=650,
                              plot_height=450)

    p.xgrid.grid_line_color = 'white'
    p.ygrid.grid_line_color = 'white'
    p.xaxis.axis_label = xlabel
    p.yaxis.axis_label = ylabel

    # Add the hover tool
    p.add_tools(hover)
    p, source1 = add_points(p, to_plot_not, 'log2FoldChange',
                            'pvalue', color='#1a9641', maxvalue=maxvalue)
    p, source2 = add_points(p, to_plot_yes, 'log2FoldChange', 'pvalue',
                            color='#fc8d59', alpha=0.6, outline=True, maxvalue=maxvalue)
    if showlabels:
        labels = LabelSet(x='log2FoldChange', y='transformed_q', text_font_size='7pt', text="gene_id", level="glyph",
                          x_offset=5, y_offset=5, source=source2, render_mode='canvas')
        p.add_layout(labels)
    if searchbox:
        text = TextInput(title="text", value="gene")
        text.js_on_change('value', CustomJS(
            args=dict(source=source1), code="""
      var data = source.data
      var value = cb_obj.value
      var gene_id = data.gene_id
      var a = -1
      for (i=0; i < gene_id.length; i++) {
          if ( gene_id[i]===value ) { a=i; console.log(i); data.size[i]=7; data.alpha[i]=1; data.color[i]='#fc8d59' }
      }
      source.data = data
      console.log(source)
      console.log(cb_obj)
      source.change.emit()
      console.log(source)
      """))
        p = column(text, p)
    if folder:
      save(p, folder + title.replace(' ', "_") + "_volcano.html")
      #export_png(p, filename=folder + title.replace(' ', "_") + "_volcano.png")
    try:
        show(p)
    except:
        show(p)
    return p


def add_points(p, df1, x, y, color='blue', alpha=0.2, outline=False, maxvalue=100):
    """parts of volcano plot"""
    # Define colors in a dictionary to access them with
    # the key from the pandas groupby funciton.
    df = df1.copy()
    transformed_q = -df[y].apply(np.log10).values
    transformed_q[transformed_q == np.inf] = maxvalue
    df['transformed_q'] = transformed_q
    df['color'] = color
    df['alpha'] = alpha
    df['size'] = 7
    source1 = bokeh.models.ColumnDataSource(df)

    # Specify data source
    p.scatter(x=x, y='transformed_q', size='size',
              alpha='alpha', source=source1,
              color='color', name='circles')
    if outline:
        p.scatter(x=x, y='transformed_q', size=7,
                  alpha=1,
                  source=source1, color='black',
                  fill_color=None, name='outlines')

    # prettify
    p.background_fill_color = "#DFDFE5"
    p.background_fill_alpha = 0.5
    return p, source1


def selector(df, valtoextract=[], logfoldtohighlight=0.15, pvaltohighlight=0.1, minlogfold=0.15, minpval=0.1):
    """Part of Volcano plot: A function to separate tfs from everything else"""
    toshow = (df.pvalue < minpval) & (abs(df.log2FoldChange) > minlogfold)
    df = df[toshow]
    sig = (df.pvalue < pvaltohighlight) & (
        abs(df.log2FoldChange) > logfoldtohighlight)
    if valtoextract:
        not_tf = (~df.gene_id.isin(valtoextract))
        is_tf = (df.gene_id.isin(valtoextract))
        to_plot_not = df[~sig | not_tf]
        to_plot_yes = df[sig & is_tf]
    else:
        to_plot_not = df[~sig]
        to_plot_yes = df[sig]
    return to_plot_not, to_plot_yes

# What pops up on hover?


def correlationMatrix(data, names, colors=None, pvals=None, maxokpval=10**-9, other=None, title="correlation Matrix", dataIsCorr=False,
                          invert=False, size=40, folder='', interactive=False, maxval=None, minval=None):
    """
    Make an interactive correlation matrix from an array using bokeh

    Args:
    -----
      data: arrayLike of int / float/ bool of size(names*val) or (names*names)
      names: list[str] of names for each rows
      colors: list[int] of size(names) a color for each names (good to display clusters)
      pvals: arraylike of int / float/ bool of size(names*val) or (names*names) with the corresponding pvalues
      maxokpval: float threshold when pvalue is considered good. otherwise lowers the size of the square
        until 10**-3 when it disappears
      other: arrayLike of int / float/ bool of size(names*val) or (names*names), an additional information
        matrix that you want ot display with opacity whereas correlations willl be displayed with
      title: str the plot title
      dataIsCorr: bool if not true, we will compute the corrcoef of the data array
      invert: bool whether or not to invert the matrix before running corrcoef
      size: int the plot size
      folder: str of folder location where to save the plot, won't save if empty
      interactive: bool whether or not to make the plot interactive (else will use matplotlib)
      maxval: float clamping coloring up to maxval
      minval: float clamping coloring down to minval

    Returns:
    -------
      the bokeh object if interactive else None

    """
    if not dataIsCorr:
        print("computing correlations")
        data = np.corrcoef(np.array(data) if not invert else np.array(data).T)
    else:
        data = np.array(data)
    regdata = data.copy()
    if minval is not None:
        data[data<minval]=minval
    if maxval is not None:
        data[data > maxval] = maxval
    data=data/data.max()
    TOOLS = "hover,crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,save"
    xname = []
    yname = []
    color = []
    alpha = []
    height = []
    width = []
    if type(colors) is list:
        print('we are assuming you want to display clusters with colors')
    elif other is not None:
        print('we are assuming you want to display the other of your correlation with opacity')
    if pvals is not None:
        print('we are assuming you want to display the pvals of your correlation with size')
        regpvals = pvals.copy()
        u = pvals<maxokpval
        pvals[~u] = np.log10(1/pvals[~u])
        pvals = pvals/pvals.max()
        pvals[u]=1
    if interactive:
        xname = []
        yname = []
        color = []
        for i, name1 in enumerate(names):
            for j, name2 in enumerate(names):
                xname.append(name1)
                yname.append(name2)
                if pvals is not None:
                    height.append(max(0.1, min(0.9, pvals[i, j])))
                    color.append(cc.coolwarm[int((data[i, j]*128)+127)])
                    alpha.append(min(abs(data[i, j]), 0.9))
                elif other is not None:
                    color.append(cc.coolwarm[int((data[i, j]*128)+127)])
                    alpha.append(max(min(other[i, j], 0.9),0.1) if other[i, j]!=0 else 0)
                else:
                    alpha.append(min(abs(data[i, j]), 0.9))
                if colors is not None:
                    if type(colors) is list:
                        if colors[i] == colors[j]:
                            color.append(
                                Category10[10][colors[i]])
                        else:
                            color.append('lightgrey')

                elif pvals is None and other is None:
                    color.append('grey' if data[i, j]
                                 > 0 else Category20[3][2])
        print(regdata.max())
        if pvals is not None:
            width = height.copy()
            data = dict(
                xname=xname,
                yname=yname,
                colors=color,
                alphas=alpha,
                data=regdata.ravel(),
                pvals=regpvals.ravel(),
                width=width,
                height=height
            )
        else:
            data = dict(
                xname=xname,
                yname=yname,
                colors=color,
                alphas=alpha,
                data=data.ravel()
            )
        tt = [('names', '@yname, @xname'), ('value', '@data')]
        if pvals is not None:
            tt.append(('pvals','@pvals'))
        p = figure(title=title if title is not None else "Correlation Matrix",
                   x_axis_location="above", tools=TOOLS,
                   x_range=list(reversed(names)), y_range=names,
                   tooltips=tt)

        p.plot_width = 800
        p.plot_height = 800
        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        p.axis.major_label_text_font_size = "5pt"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = np.pi / 3

        p.rect('xname', 'yname', width = 0.9 if not width else 'width',
                height = 0.9 if not height else 'height', source=data,
               color='colors', alpha='alphas', line_color=None,
               hover_line_color='black', hover_color='colors')
        save(p, folder + title.replace(' ', "_") + "_correlation.html")
        #export_png(p, filename=folder + title.replace(' ', "_") + "_correlation.png")
        try:
            show(p)
        except:
            show(p)
        return p  # show the plot
    else:
        plt.figure(figsize=(size, 200))
        plt.title('the correlation matrix')
        plt.imshow(data)
        plt.savefig(title + "_correlation.pdf")
        plt.show()


def venn(inp, names, title="venn", folder=''):
    """
    Plots a venn diagram using the pyvenn package

    Args:
    -----
      inp: list[set()] of sets of values (e.g. [(1,2,3,4),(2,3),(1,3,4,5)])
      names: list[str] of the name of each leaf
      title: str the plot title
      folder: str of location where to save the plot, won't save if empty
    """
    labels = pyvenn.get_labels(inp, fill=['number', 'logic'])
    if len(inp) == 2:
        fig, ax = pyvenn.venn2(labels, names=names)
    elif len(inp) == 3:
        fig, ax = pyvenn.venn3(labels, names=names)
    elif len(inp) == 4:
        fig, ax = pyvenn.venn4(labels, names=names)
    elif len(inp) == 5:
        fig, ax = pyvenn.venn5(labels, names=names)
    elif len(inp) == 6:
        fig, ax = pyvenn.venn6(labels, names=names)
    else:
        raise ValueError('need to be between 2 to 6')
    ax.set_title(title)
    if folder:
        fig.savefig(folder + title + '_venn.pdf')
    fig.show()
    plt.pause(0.1)


def mergeImages(images, outputpath):
    """
    will merge a set of images in python

    Args:
    -----
      images: list of image filepath
      outputpath: where to save the resulting merger
    """
    images = list(map(Image.open, images))
    widths, heights = zip(*(i.size for i in images))

    total_width = max(widths)
    max_height = sum(heights)

    new_im = Image.new('RGB', (total_width, max_height))

    y_offset = 0
    for im in images:
        new_im.paste(im, (0, y_offset))
        y_offset += im.size[1]

    new_im.save(outputpath)


def addTextToImage(image, text, outputpath, xy=(0, 0), color=(0, 0, 0), fontSize=64):
    """
    will add some text to an image in python

    Args:
    ----
      image: the image filepath
      text: the text to write
      outputpath: the location of the resulting image
      xy: the location of the text
      color: tuple(a,b,c) a tuple of 3 ints between 0 and 256
      fontSize: an int for the font size
    """
    # adds black text to the upper left by default, Arial size 64
    img = Image.open(image)
    draw = ImageDraw.Draw(img)
    # the below file path assumes you're operating macOS
    font = ImageFont.truetype("/Library/Fonts/Arial.ttf", fontSize)
    draw.text(xy, text, color, font=font)
    img.save(outputpath)


def SOMPlot(net, size, colnames, minweight=0.1, distq1=0.535, distq2=0.055, distr=0.2, folder=""):
    """
    makes an interactive plot from a SOM from the SimpSOM package

    a tool that uses simpSOM's package output (which produces self organizing maps),
    to plot its output in an interactive fashion

    Args:
    -----
      net: SIMPSOM's net object, output from somplot
      size: float the size of the plot
      distq1: float a value to adjust to considering the numbe of nodes (controls spacing beween cols of nodes)
      distq2: float a value to adjust to considering the numbe of nodes (controls spacing beween rows of nodes)
      distr: float a value to adjust to considering the numbe of nodes (controls nodes size spacing)
      folder: str of location where to save the plot, won't save if empty
    """
    diffs = net.diff_graph(show=False, returns=True)
    somnodes = {'r':[],'q':[],'c':diffs,'features':[]}
    for i, node in enumerate(net.nodeList):
        somnodes['q'].append(node.pos[0]+(i%size)*distq1+(i//size)*distq2)
        somnodes['r'].append(-node.pos[1]-(i%size)*distr)
        somnodes['features'].append([colnames[i] for i in np.argsort(
            node.weights) if abs(node.weights[i]) > minweight])
    somnodes=pd.DataFrame(somnodes)
    for i, v in somnodes.iterrows():
        tot=""
        for e, j in enumerate(v.features):
            if e%5==4:
                tot+='\n'
            tot += " "+str(j)
        somnodes.loc[i, 'features'] = tot
    #interactive SOM with features with highest importance to the nodes, displayed when hovering
    bigScatter(somnodes, precomputed=True, features=True, binsize=1, title='Cobinding SOM cluster of '+str(size), folder=folder)
