#!/usr/bin/env python3
# Author: Dominik Moritz, originally for spacegraphcats.

class GmlWriter:
    """Similar to the writer for gxt above but for gml."""

    def __init__(self, file, vertex_attributes=None, edge_attributes=None, directed=False):
        """Initialize graph writer."""
        self.file = file

        if vertex_attributes is not None:
            self.vertex_attributes = vertex_attributes
        if edge_attributes is not None:
            self.edge_attributes = edge_attributes

        if directed:
            self._write('graph [\n   directed 1\n')
        else:
            self._write('graph [\n   directed 0\n')

    def _write(self, string):
        self.file.write(string)

    def _quote(self, value):
        if isinstance(value, str):
            return '"{}"'.format(value)
        return value

    def add_vertex(self, id, size, attribute_values=[], vertex_attributes=None):
        """Add a vertex to the output."""
        if not hasattr(self, 'vertex_attributes'):
            self.vertex_attributes = vertex_attributes

        self._write('  node [\n')
        self._write('    id {}\n'.format(id))
        self._write('    size {}\n'.format(size))
        for k, v in zip(self.vertex_attributes, attribute_values):
            self._write('    {} {}\n'.format(k, self._quote(v)))
        self._write('  ]\n')

    def add_edge(self, src, dest, attribute_values=[], edge_attributes=None):
        """Add an edge to the output."""
        if not hasattr(self, 'edge_attributes'):
            self.edge_attributes = edge_attributes

        self._write('  edge [\n')
        self._write('    source {}\n'.format(src))
        self._write('    target {}\n'.format(dest))
        for k, v in zip(self.edge_attributes, attribute_values):
            self._write('    {} {}\n'.format(k, self._quote(v)))
        self._write('  ]\n')

    def done(self):
        """Call when done."""
        self._write(']\n')


class DotWriter:
    """Similar to the writer for gxt above but for dot."""

    def __init__(self, file):
        """Initialize graph writer."""
        self.file = file

        self._write('graph G {\n')

    def _write(self, string):
        self.file.write(string)

    def add_vertex(self, id):
        """Add a vertex to the output."""
        self._write('  {};\n'.format(id))

    def add_edge(self, src, dest):
        """Add an edge to the output."""
        self._write('  {} -- {};\n'.format(src, dest))

    def done(self):
        """Call when done."""
        self._write('}\n')
