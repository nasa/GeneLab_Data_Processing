__author__ = 'Jonathan Rubin'

#THIS SCRIPT IS DEPRECATED (JDR: 8/5/18). It currently spits out an html file with an interactive MA and volcano plot but the newer galaxy_mode.py does the same
#but with a much better html 'report' including a PCA plot and a heatmap
  
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['image.cmap'] = 'jet'
import os, config, math, mpld3, warnings, json, pylab
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.stats import gaussian_kde

#mpld3 hack to make it run with matplotlib v2.2.2
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        import numpy as np
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

from mpld3 import _display
_display.NumpyEncoder = NumpyEncoder

#This is a mpld3 plugin class to add a slider bar. Not currently in use (6/27/18) but potentially useful for a padj cutoff slider bar. Would need to
#make a web-app though and this might be out of scope for this project.
class SliderView(mpld3.plugins.PluginBase):
    """ Add slider and JavaScript / Python interaction. """

    JAVASCRIPT = """
    mpld3.register_plugin("sliderview", SliderViewPlugin);
    SliderViewPlugin.prototype = Object.create(mpld3.Plugin.prototype);
    SliderViewPlugin.prototype.constructor = SliderViewPlugin;
    SliderViewPlugin.prototype.requiredProps = ["idline", "callback_func"];
    SliderViewPlugin.prototype.defaultProps = {}

    function SliderViewPlugin(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    SliderViewPlugin.prototype.draw = function(){
      var line = mpld3.get_element(this.props.idline);
      var callback_func = this.props.callback_func;

      var div = d3.select("#" + this.fig.figid);

      // Create slider
      div.append("input").attr("type", "range").attr("min", 0).attr("max", 1).attr("step", 0.01).attr("value", 1)
          .on("change", function() {
              var command = callback_func + "(" + this.value + ")";
              console.log("running "+command);
              var callbacks = { 'iopub' : {'output' : handle_output}};
              var kernel = IPython.notebook.kernel;
              kernel.execute(command, callbacks, {silent:false});
          });

      function handle_output(out){
        //console.log(out);
        var res = null;
        // if output is a print statement
        if (out.msg_type == "stream"){
          res = out.content.data;
        }
        // if output is a python object
        else if(out.msg_type === "pyout"){
          res = out.content.data["text/plain"];
        }
        // if output is a python error
        else if(out.msg_type == "pyerr"){
          res = out.content.ename + ": " + out.content.evalue;
          alert(res);
        }
        // if output is something we haven't thought of
        else{
          res = "[out type not implemented]";  
        }

        // Update line data
        line.data = JSON.parse(res);
        line.elements()
          .attr("d", line.datafunc(line.data))
          .style("stroke", "black");

       }

    };
    """

    def __init__(self, line, callback_func):
        self.dict_ = {"type": "sliderview",
                      "idline": mpld3.utils.get_id(line),
                      "callback_func": callback_func}

#This is in relation to the above class. Not exactly sure how it's supposed to be implemented but I found this online.
def updateSlider(pval_cut):
    return pval_cut

#This function is the one currently in use (JDR 6/20/18,6/26/18) and uses the python package mpld3 to plot interactive MA and volcano plots.
#I found this one to be the most versatile and least laggy. One issue with this is that it REQUIRES matplotlib 1.3.1 (or at least does not work with
#matplotlib 2.2.2). The reason is that there is some issue with matplotlib 2.2.2 figures not being json serializable. This is really an mpld3 bug but
#if you just use matplotlib 1.3.1, it works fine.
def differential_visualize(rawdata_out,GLDS):
    #In this section of the code, the style of labels is defined. In this case, we are using a table with non_sig hits being blue and sig hits being red
    css = """
        table
        {
          border-collapse: collapse;
        }
        th
        {
          color: #ffffff;
          background-color: #000000;
        }
        td
        {
          background-color: #cccccc;
          color: #ffffff;
        }
        table, th, td
        {
          font-family:Arial, Helvetica, sans-serif;
          border: 1px solid white;
          text-align: right;
          opacity: 0.9;
        }
        """

    sig = """<table>
        <tr>
            <th></th>
            <th>{gene}</th>
        </tr>
        <tr>
            <td style="background-color: #000000;">x</td>
            <td style="background-color: #ff0000;">{x}</td>
        </tr>
        <tr>
            <td style="background-color: #000000;">y</td>
            <td style="background-color: #ff0000;">{y}</td>
        </tr>
        <tr>
            <td style="background-color: #000000;">p</td>
            <td style="background-color: #ff0000;">{pval}</td>
        </tr>
        </table>"""

    non_sig = """<table>
        <tr>
            <th></th>
            <th>{gene}</th>
        </tr>
        <tr>
            <td style="background-color: #000000;">x</td>
            <td style="background-color: #0000ff;">{x}</td>
        </tr>
        <tr>
            <td style="background-color: #000000;">y</td>
            <td style="background-color: #0000ff;">{y}</td>
        </tr>
        <tr>
            <td style="background-color: #000000;">p</td>
            <td style="background-color: #0000ff;">{pval}</td>
        </tr>
        </table>"""


    #In this section, the differential expression file is parsed and all necessary lists are generated (to avoid looping through multiple times)
    foldChange = list()
    averageExpression = list()
    adjustedPvalue = list()
    geneName = list()
    log10pval = list()
    scattersigx = list()
    scattersigy = list()
    volcanosigx = list()
    volcanosigy = list()
    scatterlabels = list()
    scattersiglabels = list()
    volcanolabels = list()
    volcanosiglabels = list()
    cell_text = list()
    condition1,condition2,pval_cut = config.visualize.split(',')
    pval_cut = float(pval_cut)
    diffExp_file = os.path.join(rawdata_out,GLDS + "_microarray_DGE.txt")
    with open(diffExp_file) as F:
        header = F.readline().strip('\n').split('\t')
        fc_index = [i for i in range(len(header)) if 'FC' in header[i]][0]+1
        exp_index = [i for i in range(len(header)) if 'Exp' in header[i]][0]+1
        p_index = [i for i in range(len(header)) if 'adj' in header[i]][0]+1
        for line in F:
            linelist = line.strip('\n').split('\t')
            gene = linelist[0].strip('"')
            fc = float(linelist[fc_index])
            exp = float(linelist[exp_index])
            pval = float(linelist[p_index])
            fc_short = "%.3f" % fc
            exp_short = "%.3f" % exp
            pval_short = "%.3f" % pval
            geneName.append(gene)
            foldChange.append(fc)
            averageExpression.append(exp)
            adjustedPvalue.append(pval)
            if pval < pval_cut:
                cell_text.append([gene,exp,fc,pval])
                scattersigx.append(exp)
                scattersigy.append(fc)
                scattersiglabels.append(sig.format(gene=gene,x=exp_short,y=fc_short,pval=pval_short))
                try:
                    l10p = -math.log(pval,10)
                    log10pval.append(l10p)
                    volcanosigy.append(l10p)
                    volcanosigx.append(fc)
                    volcanosiglabels.append(sig.format(gene=gene,x=fc_short,y="%.3f" % -math.log(pval,10),pval=pval_short))
                except ValueError:
                    print "Error: Zero adjusted p-value encountered, cannot display in volcano plot.."
                    
            else:
                log10pval.append(-math.log(pval,10))
                scatterlabels.append(non_sig.format(gene=gene,x=exp_short,y=fc_short,pval="%.3f" % pval))
                volcanolabels.append(non_sig.format(gene=gene,x=fc_short,y="%.3f" % -math.log(pval,10),pval=pval_short))

            
    #In this section the matplotlib figure is initialized and the MA-plot is created
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        F = plt.figure(figsize=(14,6.5))
        gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1])
        ax0 = F.add_subplot(gs[0])
        ax0.grid(color='black', linestyle='dashed')
        x = averageExpression
        y = foldChange
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        idx = np.argsort(z)
        x, y, z = [x[i] for i in idx], [y[i] for i in idx], [z[i] for i in idx]
        scatter = ax0.scatter(x=x,y=y,c=z,s=100,edgecolor="")
        sigscatter = ax0.scatter(scattersigx,scattersigy,c='r',s=100,edgecolor="")
        ax0.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
        ax0.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
        ax0.set_title("Microarray MA Plot", size=25)
        ax0.set_ylabel("Log2 Fold Change ("+condition1+" - "+condition2+")", size=18)
        ax0.set_xlabel("Average Expression", size=18)


        #Here we create the Volcano plot
        ax1 = F.add_subplot(gs[1])
        x = foldChange
        y = log10pval
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        idx = np.argsort(z)
        x, y, z = [x[i] for i in idx], [y[i] for i in idx], [z[i] for i in idx]
        volcano = ax1.scatter(x=x,y=y,c=z,s=100,edgecolor="")
        sigvolcano = ax1.scatter(volcanosigx,volcanosigy,c='r',s=100,edgecolor="")
        ax1.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
        ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
        ax1.set_title("Microarray Volcano Plot", size=25)
        ax1.set_ylabel("-Log10 P-value", size=18)
        ax1.set_xlabel("Log2 Fold Change ("+condition1+" - "+condition2+")", size=18)
        ax1.grid(color='black', linestyle='dashed')
        ax1.set_ylim(bottom=0)

        #This line adjusts the whitespace around the subplots since I noticed there was a lot of wasted space
        F.subplots_adjust(left=0.05,right=0.95,hspace = 0.05, wspace = 0.15)

        #This is the bulk of the mpld3 code, basically we're creating 'tooltips' which allows us to have the interactive labels
        tooltip = mpld3.plugins.PointHTMLTooltip(scatter, labels=scatterlabels[:int(0.5*len(scatterlabels))], css=css)
        tooltip2 = mpld3.plugins.PointHTMLTooltip(sigscatter, labels=scattersiglabels, css=css)
        tooltip3 = mpld3.plugins.PointHTMLTooltip(volcano, labels=volcanolabels[:int(0.5*len(volcanolabels))], css=css)
        tooltip4 = mpld3.plugins.PointHTMLTooltip(sigvolcano, labels=volcanosiglabels, css=css)
        # mpld3.plugins.connect(F, tooltip, tooltip2, tooltip3, tooltip4, SliderView(scatter, callback_func="updateSlider"))

        #This connects our plots together
        mpld3.plugins.connect(F, tooltip, tooltip2, tooltip3, tooltip4)

        #Here we save both a png version of the plot (non-interactive) and the interactive html version of the plot
        plt.savefig(os.path.join(rawdata_out,GLDS+'_visualization_results.png'))
        mpld3.save_html(F,os.path.join(rawdata_out,GLDS+'_visualization_results.html'))
        plt.close(F)

    #This section of the code creates an html table of significant genes
    with open(os.path.join(rawdata_out,GLDS+'_significant_genes_'+str(pval_cut)+'.html'),'w') as sigGenes_file:
        sigGenes_file.write("""<!DOCTYPE html>
            <html>
            <head>
            <title>List of Significant Genes</title>
            <style>
            table {
                font-family: arial, sans-serif;
                border-collapse: collapse;
                width: 100%;
            }

            td, th {
                border: 1px solid #dddddd;
                text-align: left;
                padding: 8px;
            }

            tr:nth-child(even) {
                background-color: #dddddd;
            }
            </style>
            </head>
            <body style="width: 1300px; overflow:scroll">
                <h1>List of Significant Genes</h1>
            <div>
                <div style="float: middle; width: 1300px; overflow:scroll; padding-bottom:25px; padding-top:25px">
                    <table> 
                        <tr>
                            <th>Gene</th>
                            <th>Average Expression</th> 
                            <th>Log Fold Change</th>
                            <th>Adjusted P-value</th>
                        </tr>""")
        for row in cell_text:
            name,exp,fc,pval = row
            sigGenes_file.write("""
                        <tr>
                            <td>"""+name+"""</td>
                            <td>"""+str(exp)+"""</td>
                            <td>"""+str(fc)+"""</td>
                            <td>"""+str(pval)+"""</td>
                        </tr>""")
        sigGenes_file.write("""        </table>
                </div>
            </div>
            
            </body>
            </html>""")


    #Finally this short section appends a link to the bottom of the graph html that will go directly to the list of significant genes
    with open(os.path.join(rawdata_out,GLDS+'_visualization_results.html'),'a') as html_file:
        html_file.write('<b>There were <a style="font-size: 20" href="'+GLDS+'_significant_genes_'+str(pval_cut)+'.html">'+str(len(scattersigx))+' significant genes</a> called with p-adj < '+str(pval_cut)+'</b>')


#This is not currently in use (6/26/18: JDR) it is a function to create a plotly figure which looks fairly nice but is pretty laggy
def MA_plotly(diffExp_file):
    import plotly.plotly as py
    import plotly
    import plotly.graph_objs as go  
    sig = list()
    sig_text = list()
    non_sig = list()
    non_sig_text = list()
    pval_cut = 0.1
    with open(diffExp_file) as F:
        header = F.readline().strip('\n').split('\t')
        fc_index = [i for i in range(len(header)) if 'FC' in header[i]][0]+1
        exp_index = [i for i in range(len(header)) if 'Exp' in header[i]][0]+1
        p_index = [i for i in range(len(header)) if 'adj' in header[i]][0]+1
        for line in F:
            linelist = line.strip('\n').split('\t')
            FC = linelist[fc_index]
            Exp = linelist[exp_index]
            pval = float(linelist[p_index])
            if pval < pval_cut:
                sig.append((Exp,FC))
                sig_text.append(linelist[0].strip('"'))
            else:
                non_sig.append((Exp,FC))
                non_sig_text.append(linelist[0].strip('"'))


    trace_comp0 = go.Scatter(
        x=[x for x,y in non_sig],
        y=[y for x,y in non_sig],
        mode='markers',
        marker=dict(size=12,
                    color="navy"
                   ),
        name='Non-Significant Probes',
        text=non_sig_text,
        )

    trace_comp1 = go.Scatter(
        x=[x for x,y in sig],
        y=[y for x,y in sig],
        mode='markers',
        marker=dict(size=12,
                    color="red"
                   ),
        name='Significant Probes',
        text=sig_text,
            )

    data_comp = [trace_comp0, trace_comp1]
    layout_comp = go.Layout(
        title='Microarray MA-Plot',
        hovermode='closest',
        xaxis=dict(
            title='Average Expression',
            ticklen=5,
            zeroline=False,
            gridwidth=2,
        ),
        yaxis=dict(
            title='Log Fold Change',
            ticklen=5,
            gridwidth=2,
            zeroline=True
        ),
    )
    fig_comp = go.Figure(data=data_comp, layout=layout_comp)
    plotly.offline.plot(fig_comp, filename='/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/MA-Plot.html',
        auto_open=False)

#This is not currently in use (6/26/18: JDR) - it is a function to create pygal interactive figures but it's also kind of laggy and not great looking.
def MA_pygal(diffExp_file):
    import pygal
    from pygal import config
    foldChange = list()
    averageExpression = list()
    adjustedPvalue = list()
    geneName = list()
    pval_cut = 0.1
    with open(diffExp_file) as F:
        header = F.readline().strip('\n').split('\t')
        fc_index = [i for i in range(len(header)) if 'FC' in header[i]][0]+1
        exp_index = [i for i in range(len(header)) if 'Exp' in header[i]][0]+1
        p_index = [i for i in range(len(header)) if 'adj' in header[i]][0]+1
        for line in F:
            linelist = line.strip('\n').split('\t')
            geneName.append(linelist[0].strip('"'))
            foldChange.append(float(linelist[fc_index]))
            averageExpression.append(float(linelist[exp_index]))
            adjustedPvalue.append(float(linelist[p_index]))

    xy_chart = pygal.XY(stroke=False,title="Microarray MA-Plot",x_title="Log Fold Change",y_title="Average Expression",show_legend=False)
    xy_chart.add('Sig',[{'value':(x,y),'label':g} for g,x,y,p in zip(geneName,averageExpression,foldChange,adjustedPvalue) if p < pval_cut])
    xy_chart.add('Non-Sig',[{'value':(x,y),'label':g} for g,x,y,p in zip(geneName,averageExpression,foldChange,adjustedPvalue) if p >= pval_cut])
    xy_chart.render_to_file('/Users/jonathanrubin/Google Drive/NASA/home/batch_out/GLDS-4/microarray/MA-Plot_pygal.svg')

