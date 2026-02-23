import networkx as nx
import matplotlib.pyplot as plt
import pyvis.network as pv


def draw_graph3(networkx_graph,notebook=True,output_filename='graph.html',show_buttons=True,only_physics_buttons=False,
                height=None,width=None,bgcolor=None,font_color=None,pyvis_options=None):
    """
    This function accepts a networkx graph object,
    converts it to a pyvis network object preserving its node and edge attributes,
    and both returns and saves a dynamic network visualization.
    Valid node attributes include:
        "size", "value", "title", "x", "y", "label", "color".
        (For more info: https://pyvis.readthedocs.io/en/latest/documentation.html#pyvis.network.Network.add_node)
    Valid edge attributes include:
        "arrowStrikethrough", "hidden", "physics", "title", "value", "width"
        (For more info: https://pyvis.readthedocs.io/en/latest/documentation.html#pyvis.network.Network.add_edge)
    Args:
        networkx_graph: The graph to convert and display
        notebook: Display in Jupyter?
        output_filename: Where to save the converted network
        show_buttons: Show buttons in saved version of network?
        only_physics_buttons: Show only buttons controlling physics of network?
        height: height in px or %, e.g, "750px" or "100%
        width: width in px or %, e.g, "750px" or "100%
        bgcolor: background color, e.g., "black" or "#222222"
        font_color: font color,  e.g., "black" or "#222222"
        pyvis_options: provide pyvis-specific options (https://pyvis.readthedocs.io/en/latest/documentation.html#pyvis.options.Options.set)
    """

    # import
    from pyvis import network as net

    # make a pyvis network
    network_class_parameters = {"notebook": notebook, "height": height, "width": width, "bgcolor": bgcolor, "font_color": font_color}
    pyvis_graph = net.Network(**{parameter_name: parameter_value for parameter_name, parameter_value in network_class_parameters.items() if parameter_value})

    # for each node and its attributes in the networkx graph
    for node,node_attrs in networkx_graph.nodes(data=True):
        pyvis_graph.add_node(node,**node_attrs)

    # for each edge and its attributes in the networkx graph
    for source,target,edge_attrs in networkx_graph.edges(data=True):
        # if value/width not specified directly, and weight is specified, set 'value' to 'weight'
        if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
            # place at key 'value' the weight of the edge
            edge_attrs['value']=edge_attrs['weight']
        # add the edge
        pyvis_graph.add_edge(source,target,**edge_attrs)

    # turn buttons on
    if show_buttons:
        if only_physics_buttons:
            pyvis_graph.show_buttons(filter_=['physics'])
        else:
            pyvis_graph.show_buttons()

    # pyvis-specific options
    if pyvis_options:
        pyvis_graph.set_options(pyvis_options)

    # return and also save
    pyvis_graph.force_atlas_2based()
    pyvis_graph.set_edge_smooth('dynamic')
    #pyvis_graph.set_options()


    return pyvis_graph.show(output_filename)


G2 = nx.Graph()

from Bio.SubsMat import MatrixInfo

G2.add_node('G', shape='circularImage', image='glycine.png', color='grey', hover='red')
G2.add_node('A', shape='circularImage', image='alanine.png', color='grey')
G2.add_node('V', shape='circularImage', image='valine.png', color='grey')
G2.add_node('I', shape='circularImage', image='isoleucine.png', color='grey')
G2.add_node('L', shape='circularImage', image='leucine.png', color='grey')
G2.add_node('P', shape='circularImage', image='proline.png', color='grey')
G2.add_node('F', shape='circularImage', image='phenylalanine.png', color='grey')
G2.add_node('W', shape='circularImage', image='tryptophan.png', color='grey')
G2.add_node('Y', shape='circularImage', image='tyrosine.png', color='grey')
G2.add_node('D', shape='circularImage', image='aspartic_acid.png', color='grey')
G2.add_node('E', shape='circularImage', image='glutamic_acid.png', color='grey')
G2.add_node('R', shape='circularImage', image='arginine.png', color='grey')
G2.add_node('H', shape='circularImage', image='histidine.png', color='grey')
G2.add_node('K', shape='circularImage', image='lysine.png', color='grey')
G2.add_node('S', shape='circularImage', image='serine.png', color='grey')
G2.add_node('T', shape='circularImage', image='threonine.png', color='grey')
G2.add_node('C', shape='circularImage', image='cysteine.png', color='grey')
G2.add_node('M', shape='circularImage', image='methionine.png', color='grey')
G2.add_node('N', shape='circularImage', image='asparagine.png', color='grey')
G2.add_node('Q', shape='circularImage', image='glutamine.png', color='grey')

#TODO: normalizacja score
#TODO: zmiana kolorku aktywnego node
#TODO: treshold score jako parametr
#TODO: PSSM dla pLGIC

for pair in MatrixInfo.blosum62.items():

    if not any(code in pair[0] for code in ['X', 'B', 'Z']):

        if pair[1] == -4:
            print(pair[1], pair[1]+3)
            G2.add_edges_from([(pair[0][0], pair[0][1], {'weight': pair[1] + 4, 'title': pair[1]})])
            #print(pair)

'''

G2.add_node('gly', shape='circularImage', image='glycine.png', color='grey')
G2.add_node('ala', shape='circularImage', image='alanine.png', color='grey')
G2.add_node('val', shape='circularImage', image='valine.png', color='grey')
G2.add_node('ile', shape='circularImage', image='isoleucine.png', color='grey')
G2.add_node('leu', shape='circularImage', image='leucine.png', color='grey')
G2.add_node('pro', shape='circularImage', image='proline.png', color='grey')

G2.add_node('phe', shape='circularImage', image='phenylalanine.png', color='grey')
G2.add_node('trp', shape='circularImage', image='tryptophan.png', color='grey')
G2.add_node('tyr', shape='circularImage', image='tyrosine.png', color='grey')
G2.add_node('asp', shape='circularImage', image='aspartic_acid.png', color='grey')
G2.add_node('glu', shape='circularImage', image='glutamic_acid.png', color='grey')
G2.add_node('arg', shape='circularImage', image='arginine.png', color='grey')
G2.add_node('his', shape='circularImage', image='histidine.png', color='grey')
G2.add_node('lys', shape='circularImage', image='lysine.png', color='grey')
G2.add_node('ser', shape='circularImage', image='serine.png', color='grey')
G2.add_node('thr', shape='circularImage', image='threonine.png', color='grey')
G2.add_node('cys', shape='circularImage', image='cysteine.png', color='grey')
G2.add_node('met', shape='circularImage', image='methionine.png', color='grey')
G2.add_node('asn', shape='circularImage', image='asparagine.png', color='grey')
G2.add_node('gln', shape='circularImage', image='glutamine.png', color='grey')


G2.add_edges_from([('gly', 'ala', {'title': '+ CH3', 'from': 'gly' })])
G2.add_edges_from([('gly', 'pro', {'title': '+ pyrrolidine', 'from': 'gly' })])


G2.add_edges_from([('ala', 'val', {'title': '+ branch 2xCH3', 'from': 'ala'})])
G2.add_edges_from([('ala', 'phe', {'title': '+ phenyl ring', 'from': 'ala'})])
G2.add_edges_from([('ala', 'trp', {'title': '+ indol ring', 'from': 'ala'})])
G2.add_edges_from([('ala', 'asp', {'title': '+ COOH', 'from': 'ala'})])
G2.add_edges_from([('ala', 'arg', {'title': '+ carboxyamid', 'from': 'ala'})])
G2.add_edges_from([('ala', 'his', {'title': '+ imidazole', 'from': 'ala'})])
G2.add_edges_from([('ala', 'ser', {'title': '+ OH', 'from': 'ala'})])
G2.add_edges_from([('ala', 'cys', {'title': '+ hydrosulphide', 'from': 'ala'})])
G2.add_edges_from([('ala', 'met', {'title': '+ thioether', 'from': 'ala'})])
G2.add_edges_from([('ala', 'asn', {'title': '+ carboxyamid', 'from': 'ala'})])

G2.add_edges_from([('val', 'ile', {'title': '+ CH3', 'from': 'val'})])
G2.add_edges_from([('val', 'leu', {'title': '+ CH2', 'from': 'val'})])
G2.add_edges_from([('ile', 'leu', {'title': 'move CH3', 'from': 'ile'})])


G2.add_edges_from([('phe', 'tyr', {'title': '+ OH', 'from': 'phe'})])
G2.add_edges_from([('phe', 'trp', {'title': 'change phenyl to indole', 'from': 'phe'})])
G2.add_edges_from([('phe', 'his', {'title': 'change phenyl to imidazole', 'from': 'phe'})])
G2.add_edges_from([('trp', 'his', {'title': 'change indole to imidazole', 'from': 'trp'})])


G2.add_edges_from([('asn', 'gln', {'title': '+ CH2', 'from': 'asn'})])

G2.add_edges_from([('asp', 'glu', {'title': '+ CH2', 'from': 'asp'})])

G2.add_edges_from([('ser', 'thr', {'title': '+ CH3', 'from': 'ser'})])
G2.add_edges_from([('ser', 'cys', {'title': 'change OH to hydrosulphide', 'from': 'ser'})])


'''

#draw_graph3(G2, show_buttons=True, height=800, width=1000, pyvis_options='"nodes": {    "color": {      "hover": {        "border": "rgba(231,44,233,1)"      }    }}')
draw_graph3(G2, show_buttons=False, height=800, width=1000)