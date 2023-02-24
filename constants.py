import seaborn as sns 


palette = sns.color_palette(n_colors=5)

ColorsDict = {
    'x-HSIC':palette[0], 
    'HSIC-perm':palette[1], 
    'x-dCov':palette[2], 
    'dCov-perm':palette[3]
}

markersDict = {
    'x-HSIC':'o', 
    'HSIC-perm':'^', 
    'x-dCov':'o', 
    'dCov-perm':'^'
}