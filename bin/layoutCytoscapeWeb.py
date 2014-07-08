#!/usr/bin/env python
"""layoutCytoscapeWeb.py

Author: Sam Ng
Last Updated: 2014-06-13
"""
import math, os, re, shutil, sys
from copy import deepcopy

import pandas
import pandas
from PATHMARK import Pathway

from optparse import OptionParser

html_link = """<a href="%s">%s</a>"""

html_head = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">
<html>
    
    <head>
        <title>CytoscapeWeb - %s</title>
        
        <script type="text/javascript" src="../js/min/json2.min.js"></script>
        <script type="text/javascript" src="../js/min/AC_OETags.min.js"></script>
        <script type="text/javascript" src="../js/min/cytoscapeweb.min.js"></script>
        
        <script type="text/javascript">
            window.onload = function() {
                // id of Cytoscape Web container div
                var div_id = "cytoscapeweb";
                
                // NOTE: - the attributes on nodes and edges
                //       - it also has directed edges, which will automatically display edge arrows
                var xml = '\\
"""

html_tail = """                ';
                
                // visual style we will use
                var visual_style = {
                    global: {
                        backgroundColor: "#FFFFFF"
                    },
                    nodes: {
                        shape: {
                            discreteMapper: {
                                attrName: "type",
                                entries: [
                                    { attrValue: "protein", value: "CIRCLE" },
                                    { attrValue: "complex", value: "HEXAGON" },
                                    { attrValue: "abstract", value: "ROUNDRECT" },
                                    { attrValue: "drug", value: "TRIANGLE" },
                                    { attrValue: "plot", value: "CIRCLE" },
                                    { attrValue: "legend", value: "ROUNDRECT" }
                                ]
                            }
                        },
                        borderWidth: 2,
                        borderColor: "#0000000",
                        size: {
                            passthroughMapper: { attrName: "size" }
                        },
                        color: {
                            passthroughMapper: { attrName: "color" }
                        },
                        image: {
                            passthroughMapper: { attrName: "image" }
                        },
                        labelHorizontalAnchor: "center"
                    },
                    edges: {
                        width: 2,
                        color: "#000000",
                        targetArrowShape: {
                            discreteMapper: {
                                attrName: "interaction",
                                    entries: [
                                    { attrValue: "-t|", value: "T" },
                                    { attrValue: "-t>", value: "ARROW" },
                                    { attrValue: "-a|", value: "T" },
                                    { attrValue: "-a>", value: "ARROW" },
                                    { attrValue: "-ap|", value: "T" },
                                    { attrValue: "-ap>", value: "ARROW" },
                                    { attrValue: "component>", value: "NONE" },
                                    { attrValue: "-disconnected-", value: "NONE" }
                                ]
                            }
                        },
                        style: {
                            discreteMapper: {
                                attrName: "interaction",
                                entries: [
                                    { attrValue: "-t|", value: "SOLID" },
                                    { attrValue: "-t>", value: "SOLID" },
                                    { attrValue: "-a|", value: "LONG_DASH" },
                                    { attrValue: "-a>", value: "LONG_DASH" },
                                    { attrValue: "-ap|", value: "LONG_DASH" },
                                    { attrValue: "-ap>", value: "LONG_DASH" },
                                    { attrValue: "component>", value: "LONG_DASH" },
                                    { attrValue: "-disconnected-", value: "LONG_DASH" }
                                ]
                            }
                        }
                    }
                };
                
                // initialization options
                var options = {
                    swfPath: "../swf/CytoscapeWeb",
                    flashInstallerPath: "../swf/playerProductInstall"
                };
                
                var vis = new org.cytoscapeweb.Visualization(div_id, options);
                
                vis.ready(function() {
                    // add a listener for when nodes and edges are clicked
                    vis.addListener("click", "nodes", function(event) {
                        handle_click(event);
                    })
                    .addListener("click", "edges", function(event) {
                        handle_click(event);
                    });
                    
                    function handle_click(event) {
                         var target = event.target;
                         
                         clear();
                         print("event.group = " + event.group);
                         for (var i in target.data) {
                            var variable_name = i;
                            var variable_value = target.data[i];
                            print( "event.target.data." + variable_name + " = " + variable_value );
                         }
                    }
                    
                    function clear() {
                        document.getElementById("note").innerHTML = "";
                    }
                    
                    function print(msg) {
                        document.getElementById("note").innerHTML += "<p>" + msg + "</p>";
                    }
                });

                var draw_options = {
                    // your data goes here
                    network: xml,
                    // hide edge labels
                    edgeLabelsVisible: false,
                    // let's try another layout
                    layout: "ForceDirected",
                    // set the style at initialisation
                    visualStyle: visual_style,
                    // show pan zoom
                    panZoomControlVisible: true 
                };
                
                vis.draw(draw_options);
            };
        </script>
        
        <style type= >
            * { margin: 0; padding: 0; font-family: Helvetica, Arial, Verdana, sans-serif; }
            html, body { height: 100%; width: 100%; padding: 0; margin: 0; background-color: #f0f0f0; }
            body { line-height: 1.5; color: #000000; font-size: 14px; }
            /* The Cytoscape Web container must have its dimensions set. */
            #cytoscapeweb { width: 100%; height: 80%; }
            #note { width: 100%; text-align: center; padding-top: 1em; }
            .link { text-decoration: underline; color: #0B94B1; cursor: pointer; }
        </style>
    </head>
    
    <body>
        <div id="cytoscapeweb">
            Cytoscape Web will replace the contents of this div with your graph.
        </div>
        <div id="note">
            Click node or edge for details
        </div>
    </body>
    
</html>
"""

class rgb:
    def __init__(self,r,g,b):
        self.r = int(round(r))
        self.g = int(round(g))
        self.b = int(round(b))
        
        if self.r > 255:
            self.r = 255
        elif self.r < 0:
            self.r = 0
        if self.g > 255:
            self.g = 255
        elif self.g < 0:
            self.g = 0
        if self.b > 255:
            self.b = 255
        elif self.b < 0:
            self.b = 0
    def tohex(self):
        r = self.r
        g = self.g
        b = self.b
        hexchars = '0123456789ABCDEF'
        return('#' + hexchars[r/16]+hexchars[r%16]+hexchars[g/16]+hexchars[g%16]+hexchars[b/16]+hexchars[b%16])

def getColor(value, min_value = -10, max_value = 10, min_color = rgb(0, 0, 255), max_color = rgb(255, 0, 0), zero_color = rgb(255, 255, 255)):
    try:
        fvalue = float(value)
        if fvalue != fvalue:
            raise ValueError
    except ValueError:
        color = rgb(200, 200, 200)
        return(color.tohex())
    if fvalue < 0.0:
        if fvalue < min_value:
            fvalue = 1.0
        else:
            fvalue = fvalue/min_value
        color = min_color
    else:
        if fvalue > max_value:
            fvalue = 1.0
        else:
            fvalue = fvalue/max_value
        color = max_color
    r = fvalue*float(color.r - zero_color.r) + zero_color.r
    g = fvalue*float(color.g - zero_color.g) + zero_color.g
    b = fvalue*float(color.b - zero_color.b) + zero_color.b
    color = rgb(r, g, b)
    return(color.tohex())

def getSize(value, min_value = -10, max_value = 10, min_size = 20, max_size = 100):
    try:
        fvalue = float(value)
        if fvalue != fvalue:
            raise ValueError
    except ValueError:
        return(min_size)
    if fvalue < 0.0:
        if fvalue < min_value:
            size = max_size
        else:
            size = min_size + (max_size - min_size)*(fvalue/min_value)
    else:
        if fvalue > max_value:
            size = max_size
        else:
            size = min_size + (max_size - min_size)*(fvalue/max_value)
    return(size)

def layoutCytoscapeWeb():
    ## parse arguments
    parser = OptionParser(usage = '%prog [options]')
    parser.add_option('-s', '--signature', dest='signature_file', default=None)
    parser.add_option('-p', '--pathway', dest='pathway_file', default=None)
    parser.add_option('-i', '--image', dest='image_directory', default=None)
    parser.add_option('-l', '--legend', dest='legend_file', default=None)
    parser.add_option('-o', '--output', dest='output_directory', default=None)
    options, args = parser.parse_args()
    
    assert(len(args) == 0)
    
    ## read pathway_file
    cytoscape_pathway = Pathway(options.pathway_file)
    
    ## read signature_file
    signature_frame = pandas.read_csv(options.signature_file, sep = '\t', index_col = 0)
    signatures = list(signature_frame.columns)
    for signature in signatures:
        signature_map = signature_frame[signature]
        signature_max = max(max(signature_map), 0)
        signature_min = min(min(signature_map), 0)
        ## create graphml structure
        graphml_content = """<graphml>\\
                        <key id="name" for="node" attr.name="name" attr.type="string"/>\\
                        <key id="label" for="node" attr.name="label" attr.type="string"/>\\
                        <key id="type" for="node" attr.name="type" attr.type="string"/>\\
                        <key id="color" for="node" attr.name="color" attr.type="string"/>\\
                        <key id="size" for="node" attr.name="size" attr.type="double"/>\\
                        <key id="score" for="node" attr.name="score" attr.type="double"/>\\
                        <key id="image" for="node" attr.name="image" attr.type="string"/>\\
                        <key id="interaction" for="edge" attr.name="interaction" attr.type="string"/>\\
                        <graph edgedefault="directed">\\
                        """
        node_map = {}
        for index, node in enumerate(cytoscape_pathway.nodes):
            node_map[node] = index
        
        if options.legend_file:
            graphml_content += """       <node id="%s">\\
                                <data key="name">%s</data>\\
                                <data key="label">%s</data>\\
                                <data key="type">%s</data>\\
                                <data key="color">%s</data>\\
                                <data key="size">%s</data>\\
                                <data key="score">%s</data>\\
                                <data key="image">%s</data>\\
                            </node>\\
                            """ % (len(cytoscape_pathway.nodes), 'LEGEND', '', 'legend', '#FFFFFF', 250, 0, 'img/%s' % (options.legend_file.split('/')[-1]))
        
        image_files = []
        for node in cytoscape_pathway.nodes:
            if cytoscape_pathway.nodes[node] == '__DISCONNECTED__':
                node_name = node
                node_label = ''
                node_type = 'abstract'
                node_color = '#FFFFFF'
                node_size = 15
                node_score = 0
                node_image = ''
            else:
                node_name = node
                if cytoscape_pathway.nodes[node] == 'protein':
                    node_label = node
                else:
                    node_label = ''
                node_type = cytoscape_pathway.nodes[node]
                if node in signature_map:
                    node_color = getColor(signature_map[node], min_value = signature_min, max_value = signature_max)
                    node_size = getSize(signature_map[node], min_value = signature_min, max_value = signature_max)
                    node_score = signature_map[node]
                else:
                    node_color = '#FFFFFF'
                    node_size = 20
                    node_score = 0
                node_image = ''
                if options.image_directory is not None:
                    if os.path.exists('%s/%s.png' % (options.image_directory, re.sub('[:/]', '_', node))):
                        image_files.append('%s/%s.png' % (options.image_directory, re.sub('[:/]', '_', node)))
                        node_image = 'img/%s.png' % (re.sub('[:/]', '_', node))
                        node_type = 'plot'
                        node_color = '#FFFFFF'
            
            graphml_content += """       <node id="%s">\\
                                <data key="name">%s</data>\\
                                <data key="label">%s</data>\\
                                <data key="type">%s</data>\\
                                <data key="color">%s</data>\\
                                <data key="size">%s</data>\\
                                <data key="score">%s</data>\\
                                <data key="image">%s</data>\\
                            </node>\\
                            """ % (node_map[node], node_name, node_label, node_type, node_color, node_size, node_score, node_image)
        for source in cytoscape_pathway.interactions:
            for target in cytoscape_pathway.interactions[source]:
                for interaction in cytoscape_pathway.interactions[source][target].split(';'):
                    graphml_content += """   <edge source="%s" target="%s">\\
                                    <data key="interaction">%s</data>\\
                                </edge>\\
                                """ % (node_map[source], node_map[target], interaction)
        graphml_content += """</graph>\\
                    </graphml>\\
                    """
        
        ## launch cytoscape
        if not os.path.exists(options.output_directory):
            os.mkdir(options.output_directory)
        cytoscape_directory = '%s/%s_%s' % (options.output_directory, options.pathway_file.split('/')[-1].rstrip('_pathway.tab'), signature)
        assert(not os.path.exists(cytoscape_directory))
        os.mkdir(cytoscape_directory)
        os.mkdir('%s/img' % (cytoscape_directory))
        f = open('%s/%s.html' % (cytoscape_directory, signature), 'w')
        f.write(html_head % (signature) + graphml_content + html_tail)
        f.close()
        for image_file in image_files:
            shutil.copy(image_file, '%s/img' % (cytoscape_directory))
        if options.legend_file:
            shutil.copy(options.legend_file, '%s/img' % (cytoscape_directory))
        
if __name__ == "__main__":
    layoutCytoscapeWeb()
