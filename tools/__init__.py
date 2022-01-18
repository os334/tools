#!/usr/bin/env python
# @Author: kt16
# @Date:   2020-07-30 12:52:55
# @Last Modified by:   kt16
# @Last Modified time: 2021-06-18 15:50:59
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import functools
import matplotlib
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import math
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# create a palette for umap
def cmp(palette='viridis'):
    viridis = cm.get_cmap(palette, 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    grey = np.array([215 / 256, 215 / 256, 215 / 256, 1])
    newcolors[:1, :] = grey
    newcmp = ListedColormap(newcolors)
    return(newcmp)
def get_hex(pal, n=None):
    if n is None:
        n = 5
    else:
        n = n
    cmap = cm.get_cmap(pal, n)    # PiYG
    cols = []
    for i in range(cmap.N):
        rgba = cmap(i)  # rgb2hex accepts rgb or rgba
        cols.append(matplotlib.colors.rgb2hex(rgba))
    return(cols)
def colorRampPalette(low, high, diverging = False, medium=None, n=None):
    if medium is None:
        medium_ = "#FFFFFF"
    if n is None:
        n_ = 256
    else:
        n_ = n
    def hex_to_RGB(hex):
        '''
        "#FFFFFF" -> [255,255,255]
        '''
        # Pass 16 to the integer function for change of base
        return [int(hex[i:i + 2], 16) for i in range(1, 6, 2)]
    def RGB_to_hex(RGB):
        ''' [255,255,255] -> "#FFFFFF" '''
        # Components need to be integers for hex to make sense
        if len(RGB) == 4:
            RGB = [int(x) for x in RGB[:3]]
        else:
            RGB = [int(x) for x in RGB]
        return "#" + "".join(["0{0:x}".format(v) if v < 16 else
                              "{0:x}".format(v) for v in RGB])
    def color_dict(gradient):
        '''
        Takes in a list of RGB sub-lists and returns dictionary of
        colors in RGB and hex form for use in a graphing function
        defined later on
        '''
        return {"hex": [RGB_to_hex(RGB) for RGB in gradient],
                "r": [RGB[0] for RGB in gradient],
                "g": [RGB[1] for RGB in gradient],
                "b": [RGB[2] for RGB in gradient]}
    def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
        '''
        returns a gradient list of (n) colors between
        two hex colors. start_hex and finish_hex
        should be the full six-digit color string,
        inlcuding the number sign ("#FFFFFF")
        '''
        # Starting and ending colors in RGB form
        s = hex_to_RGB(start_hex)
        f = hex_to_RGB(finish_hex)
        # Initilize a list of the output colors with the starting color
        RGB_list = [s]
        # Calcuate a color at each evenly spaced value of t from 1 to n
        for t in range(1, n):
            # Interpolate RGB vector for color at the current value of t
            curr_vector = [int(s[j] + (float(t) / (n - 1)) * (f[j] - s[j])) for j in range(3)]
            # Add it to our list of output colors
            RGB_list.append(curr_vector)
        return color_dict(RGB_list)
    if diverging:
        newcmp = ListedColormap(linear_gradient(start_hex = low, finish_hex = medium_, n = n_ / 2)['hex'] + linear_gradient(start_hex = medium_, finish_hex = high, n = n_ / 2)['hex'])
    else:
        newcmp = ListedColormap(linear_gradient(start_hex = low, finish_hex = high, n = n_)['hex'])
    return(newcmp)
# export DE results
def exportDEres(adata, key = None, column = None, filename = None, remove_mito_ribo = True):
    if key is None:
        key = 'rank_genes_groups'
    else:
        key = key
    if column is None:
        column = list(adata.uns[key]['scores'].dtype.fields.keys())[0]
    else:
        column = column
    if filename is None:
        returnDEres(adata, key= key, column = column, remove_mito_ribo = remove_mito_ribo)
    else:
        scores = pd.DataFrame(data = adata.uns[key]['scores'][column], index = adata.uns[key]['names'][column])
        lfc = pd.DataFrame(data = adata.uns[key]['logfoldchanges'][column], index = adata.uns[key]['names'][column])
        pvals = pd.DataFrame(data = adata.uns[key]['pvals'][column], index = adata.uns[key]['names'][column])
        padj = pd.DataFrame(data = adata.uns[key]['pvals_adj'][column], index = adata.uns[key]['names'][column])
        try:
            pts = pd.DataFrame(data = adata.uns[key]['pts'][column], index = adata.uns[key]['names'][column])       
        except:
            pass
        scores = scores.loc[scores.index.dropna()]
        lfc = lfc.loc[lfc.index.dropna()]
        pvals = pvals.loc[pvals.index.dropna()]
        padj = padj.loc[padj.index.dropna()]
        try:
            pts = pts.loc[pts.index.dropna()]
        except:
            pass
        try:
            dfs = [scores, lfc, pvals, padj, pts]
        except:
            dfs = [scores, lfc, pvals, padj]
        df_final = functools.reduce(lambda left, right: pd.merge(left, right, left_index = True, right_index = True), dfs)
        try:
            df_final.columns = ['scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'pts']
        except:
            df_final.columns = ['scores', 'logfoldchanges', 'pvals', 'pvals_adj']
        df_final.to_csv(filename, sep = '\t')
def returnDEres(adata, key= None, column = None, remove_mito_ribo = True):
    if key is None:
        key = 'rank_genes_groups'
    else:
        key = key
    if column is None:
        column = list(adata.uns[key]['scores'].dtype.fields.keys())[0]
    else:
        column = column
    scores = pd.DataFrame(data = adata.uns[key]['scores'][column], index = adata.uns[key]['names'][column])
    lfc = pd.DataFrame(data = adata.uns[key]['logfoldchanges'][column], index = adata.uns[key]['names'][column])
    pvals = pd.DataFrame(data = adata.uns[key]['pvals'][column], index = adata.uns[key]['names'][column])
    padj = pd.DataFrame(data = adata.uns[key]['pvals_adj'][column], index = adata.uns[key]['names'][column])
    try:
        pts = pd.DataFrame(data = adata.uns[key]['pts'][column], index = adata.uns[key]['names'][column])       
    except:
        pass
    scores = scores.loc[scores.index.dropna()]
    lfc = lfc.loc[lfc.index.dropna()]
    pvals = pvals.loc[pvals.index.dropna()]
    padj = padj.loc[padj.index.dropna()]
    try:
        pts = pts.loc[pts.index.dropna()]
    except:
        pass
    try:
        dfs = [scores, lfc, pvals, padj, pts]
    except:
        dfs = [scores, lfc, pvals, padj]
    df_final = functools.reduce(lambda left, right: pd.merge(left, right, left_index = True, right_index = True), dfs)
    try:
        df_final.columns = ['scores', 'logfoldchanges', 'pvals', 'pvals_adj', 'pts']
    except:
        df_final.columns = ['scores', 'logfoldchanges', 'pvals', 'pvals_adj']
    if remove_mito_ribo:
        df_final = df_final[~df_final.index.isin(list(df_final.filter(regex='^RPL|^RPS|^MRPS|^MRPL|^MT-', axis = 0).index))]
        df_final = df_final[~df_final.index.isin(list(df_final.filter(regex='^Rpl|^Rps|^Mrps|^Mrpl|^mt-', axis = 0).index))]
    return(df_final)
# caculate the centroid for each cluster
def calc_centroid(points):
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    centroid = (sum(x) / len(points), sum(y) / len(points))
    return(centroid)
def closest_node(node, nodes):
    nodes = np.asarray(nodes[:, :2])
    deltas = nodes - node
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    return np.argmin(dist_2)
def vmax(adata, genes, pct):
    if type(genes) is not list:
        genes = [genes]
    vm = []
    for g in genes:
        try:
            idx = adata.raw.var.index.get_loc(g)
        except:
            idx = adata.var.index.get_loc(g)
        vm.append(math.ceil(np.quantile(adata.raw.X[:, idx].toarray(), pct) * 100.0) / 100.0)
    return(vm)
def vmin(adata, genes, pct):
    if type(genes) is not list:
        genes = [genes]
    vm = []
    for g in genes:
        try:
            idx = adata.raw.var.index.get_loc(g)
        except:
            idx = adata.var.index.get_loc(g)
        vm.append(math.ceil(np.quantile(adata.raw.X[:, idx].toarray(), 1 - pct) * 100.0) / 100.0)
    return(vm)