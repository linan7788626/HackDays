import pygraphviz as PG

A = PG.AGraph(directed=True, strict=True)

A.add_edge("7th Edition", "32V")
A.add_edge("7th Edition", "Xenix")
# etc., etc.

# save the graph in dot format
A.write('ademo.dot')

# pygraphviz renders graphs in neato by default,
# so you need to specify dot as the layout engine
A.layout(prog='dot')


B = PG.AGraph("./o.dot")
