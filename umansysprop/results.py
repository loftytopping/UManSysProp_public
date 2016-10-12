# vim: set et sw=4 sts=4 fileencoding=utf-8:
#
# Copyright 2014 Dave Jones <dave@waveform.org.uk>.
#
# This file is part of umansysprop.
#
# umansysprop is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# umansysprop is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# umansysprop.  If not, see <http://www.gnu.org/licenses/>.

"""
==============================
``umansysprop.results`` Module
==============================

This module defines the classes used to encapsulate results returned by the
UManSysProp server. Each tool method on the client will return an instance of
the :class:`Result` class which in turn contains one or more :class:`Table`
instances.

Result
======

.. autoclass:: Result

Table
=====

.. autoclass:: Table
"""

from __future__ import (
    unicode_literals,
    absolute_import,
    print_function,
    division,
    )
str = type('')


import sys
import json
import itertools


class Result(list):
    """
    Represents a list of named :class:`Table` objects.

    The result of a method is represented as a sequence of tables. This class
    contains a list of :class:`Table` objects each of which may be retrieved by
    name, or by index in the list (tables with identical names are not ignored,
    but only the first table may be retrieved by name).

    .. note::

        This class has an extended string representation intended for easy
        command line debugging. Simply print an instance of the class to view a
        dump of all the tables contained within it.
    """
    def __init__(self, *tables):
        super(Result, self).__init__(tables)

    @classmethod
    def from_json(cls, obj):
        """
        This class constructor accepts a parsed JSON object (created by
        :func:`umansysprop.renderers.render_json`, or with the same structure
        produced by that function) and constructs the :class:`Result` from this
        structure.
        """
        def to_tuple(v):
            if isinstance(v, list):
                return tuple(v)
            else:
                return v

        tables = []
        for table_dict in obj:
            name = table_dict['name']
            title = table_dict['title']
            rows_title = to_tuple(table_dict['rows_title'])
            rows_unit = to_tuple(table_dict['rows_unit'])
            cols_title = to_tuple(table_dict['cols_title'])
            cols_unit = to_tuple(table_dict['cols_unit'])
            rows = []
            cols = []
            data = {}
            for datum in table_dict['data']:
                key = datum['key']
                value = datum['value']
                row_key, col_key = (to_tuple(v) for v in key)
                if row_key not in rows:
                    rows.append(row_key)
                if col_key not in cols:
                    cols.append(col_key)
                data[(row_key, col_key)] = value
            tables.append(Table(
                name, rows, cols, data=data, title=title,
                rows_title=rows_title, rows_unit=rows_unit,
                cols_title=cols_title, cols_unit=cols_unit))
        return cls(*tables)

    def __getattr__(self, name):
        for table in self:
            if table.name == name:
                return table
        return super(Result, self).__getattr__(name)

    def __str__(self):
        return '\n'.join(
            line
            for table in self
            for line in (
                table.name,
                '=' * len(table.name),
                '',
                str(table),
                '',
                )
            )


class Table(object):
    """
    Represents a single table in a :class:`Result`.

    A tool is expected to return a sequence of :class:`Table` objects in a
    :class:`Result` object. Each table has a :attr:`name` (which can be used to
    access it in the :class:`Result` object), an ordered list of keys for
    :attr:`rows` and :attr:`cols`, and a function which is used to derive the
    data for each cell. The function accepts two arguments, the row and column
    key in that order, and is expected to return a scalar value. The reason for
    constructing a table in this manner (lazy evaluation) is that it enables
    renderers to query the table structure and layout without necessarily
    calculating anything. Calculated data is cached on the assumption that such
    calculations are expensive.

    The row and column keys can be any immutable value (immutability is
    required as they will form keys in a dict at evaluation time). Keys which
    are tuples will be treated specially as renderers. For example, if each
    row key is a 2-tuple, then each row in the resulting table will have two
    row headers in two separate columns at the left of the table. This can
    aid in representing data with more than 2 dimensions in a table.

    Consider a result set keyed by values A, B, and C. The table can be
    constructed with a series of 2-tuple row keys (A, B), while the column can
    be scalar C values. The resulting table will be rendered as follows:

    +----+----+------+------+------+
    |    |    | C1   | C2   | C3   |
    +====+====+======+======+======+
    | A1 | B1 | data | data | data |
    +    +----+------+------+------+
    |    | B2 | data | data | data |
    +----+----+------+------+------+
    | A2 | B1 | data | data | data |
    +    +----+------+------+------+
    |    | B2 | data | data | data |
    +----+----+------+------+------+

    Optional attributes also exist for :attr:`title`, :attr:`rows_title`,
    :attr:`cols_title`, :attr:`rows_unit`, and :attr:`cols_unit` (these all
    default to an empty string if omitted). In the case that tuples are used
    for row or column keys, the corresponding title and unit values must be
    tuples as well.

    .. note::

        Like :class:`Result`, this class has an extended string representation
        intended for easy command line debugging. Printing an instance of
        this class will produce a human readable string representation of the
        table's row and column keys along with the calculated data.

    .. autoattribute:: as_ndarray

    .. autoattribute:: as_dataframe

    .. autoattribute:: col_dims

    .. autoattribute:: col_titles

    .. autoattribute:: cols

    .. autoattribute:: cols_iter

    .. attribute:: cols_title

        A string or tuple of strings giving the title of each column dimension.
        Note that if :attr:`col_dims` is 1, this may be either a string or
        a 1-tuple containing a string. The associated :attr:`col_titles`
        attribute may be easier to work with.

    .. attribute:: cols_unit

        A string or tuple of strings giving the units of each column dimension.
        Note that if :attr:`col_dims` is 1, this may be either a string or
        a 1-tuple containing a string. The associated :attr:`col_units`
        attribute may be easier to work with.

    .. autoattribute:: data

    .. autoattribute:: data_iter

    .. attribute:: name

        The name of the table. This is intended for scripting usage and as such
        will only ever contain a string beginning with an alphabetic character
        followed by zero or more alphanumeric characters or underscores.

    .. autoattribute:: row_dims

    .. autoattribute:: row_titles

    .. autoattribute:: rows

    .. autoattribute:: rows_iter

    .. attribute:: rows_title

        A string or tuple of strings giving the title of each row dimension.
        Note that if :attr:`row_dims` is 1, this may be either a string or
        a 1-tuple containing a string. The associated :attr:`row_titles`
        attribute may be easier to work with.

    .. attribute:: rows_unit

        A string or tuple of strings giving the units of each row dimension.
        Note that if :attr:`row_dims` is 1, this may be either a string or
        a 1-tuple containing a string. The associated :attr:`row_units`
        attribute may be easier to work with.

    .. attribute:: title

        The human readable title of the table, typically rendered in the web
        interface as the table's caption.
    """
    def __init__(
            self, name, rows, cols, func=None, data=None, title='',
            rows_title=None, cols_title=None, rows_unit=None, cols_unit=None):
        if func is None and data is None:
            raise ValueError('Either func or data must be specified')
        self._rows = tuple(rows)
        self._cols = tuple(cols)
        if not self._rows:
            raise ValueError('Table must have at least one row key')
        if not self._cols:
            raise ValueError('Table must have at least one column key')
        self._row_dims = len(self.rows[0]) if isinstance(self.rows[0], tuple) else 1
        self._col_dims = len(self.cols[0]) if isinstance(self.cols[0], tuple) else 1
        self.rows_title = self._keys_default(rows_title, self.row_dims)
        self.rows_unit = self._keys_default(rows_unit, self.row_dims)
        self.cols_title = self._keys_default(cols_title, self.col_dims)
        self.cols_unit = self._keys_default(cols_unit, self.col_dims)
        self._row_spans = self._calculate_spans(tuple(self.rows_iter), self.row_dims)
        self._col_spans = self._calculate_spans(tuple(self.cols_iter), self.col_dims)
        self._func = func
        self._data = data
        self.name = name
        self.title = title

    def _keys_default(self, value, dims):
        if value is None:
            if dims == 1:
                value = ''
            else:
                value = ('',) * dims
        if dims > 1:
            if not isinstance(value, tuple):
                raise ValueError('%r is not a tuple' % value)
            if len(value) != dims:
                raise ValueError('%r does not contain %dims elements' % (value, dims))
        return value

    def _calculate_spans(self, keys, dims):
        if not keys:
            raise ValueError('keys cannot be empty')
        if dims > 1:
            spans = []
            for dim in range(dims):
                dim_spans = []
                last_value = None
                span = 0
                for i, key in enumerate(keys):
                    value = key[dim]
                    if span == 0:
                        for j, comparison in enumerate(keys[i:]):
                            if comparison[dim] == key[dim]:
                                span += 1
                            else:
                                break
                        dim_spans.append(span)
                    else:
                        dim_spans.append(0)
                    span -= 1
                spans.append(dim_spans)
            return {key: tuple(spans[d][i] for d in range(dims)) for i, key in enumerate(keys)}
        else:
            return {key: (1,) for key in keys}

    @property
    def rows(self):
        """
        An ordered sequence of keys for the rows of the table. If
        :attr:`row_dims` is greater than one, then this is a sequence of
        tuples. These values, combined with :attr:`cols` can be used to
        index :attr:`data` in display order like so::

            for row in table.rows:
                for col in table.cols:
                    print(table.data[(row, col)])
        """
        return self._rows

    @property
    def row_dims(self):
        """
        The number of dimensions within the row keys. If this is greater than
        one, then :attr:`rows` is a sequence of tuples.
        """
        return self._row_dims

    @property
    def rows_iter(self):
        """
        Returns an iterator over :attr:`rows` where each key is returned as a
        tuple, regardless. This property is intended to make renderers simpler.
        """
        for row in self._rows:
            if self.row_dims > 1:
                yield row
            else:
                yield (row,)

    @property
    def row_titles(self):
        """
        Returns :attr:`rows_title` as a tuple, regardless. This property is
        intended to make renderers simpler.
        """
        if self.row_dims > 1:
            return self.rows_title
        else:
            return (self.rows_title,)

    @property
    def row_units(self):
        """
        Returns :attr:`rows_unit` as a tuple, regardless. This property is
        intended to make renderers simpler.
        """
        if self.row_dims > 1:
            return self.rows_unit
        else:
            return (self.rows_unit,)

    @property
    def row_spans(self):
        """
        A mapping of row keys (as tuples, as from :attr:`rows_iter`) to a tuple
        of row spans. For example, consider the following sequence of row
        keys::

            (1, 1), (1, 2), (1, 3), (2, 1), (2, 2)

        This would generate the following row spans mapping::

            {
                (1, 1): (3, 1),
                (1, 2): (0, 1),
                (1, 3): (0, 1),
                (2, 1): (2, 1),
                (2, 2): (0, 1),
                }

        Indicating that the first element of the first key should span three
        rows, and that the subsequent first elements within the spanned rows
        should not be rendered at all. This property is intended to make
        renderers that target human-readable formats simpler.
        """
        return self._row_spans

    @property
    def cols(self):
        """
        An ordered sequence of keys for the columns of the table. If
        :attr:`col_dims` is greater than one, this this is a sequence of
        tuples. These values, combined with :attr:`rows` can be used to index
        :attr:`data` in display order like so::

            for row in table.rows:
                for col in table.cols:
                    print(table.data[(row, col)])
        """
        return self._cols

    @property
    def col_dims(self):
        """
        The number of dimensions within the column keys. If this is greater
        than one, then :attr:`cols` is a sequence of tuples.
        """
        return self._col_dims

    @property
    def cols_iter(self):
        """
        Returns an iterator over :attr:`cols` where each key is returned as
        a tuple, regardless.
        """
        for col in self._cols:
            if self.col_dims > 1:
                yield col
            else:
                yield (col,)

    @property
    def col_titles(self):
        """
        Returns :attr:`cols_title` as a tuple, regardless. This property is
        intended to make renderers simpler.
        """
        if self.col_dims > 1:
            return self.cols_title
        else:
            return (self.cols_title,)

    @property
    def col_units(self):
        """
        Returns :attr:`cols_unit` as a tuple, regardless. This property is
        intended to make renderers simpler.
        """
        if self.col_dims > 1:
            return self.cols_unit
        else:
            return (self.cols_unit,)

    @property
    def col_spans(self):
        """
        A mapping of col keys (as tuples, as from :attr:`cols_iter`) to a tuple
        of col spans. See :attr:`row_spans` for an example of the mapping
        returned.
        """
        return self._col_spans

    @property
    def data(self):
        """
        The data contained within the table. This is presented as a dict
        keyed by `(row_key, col_key)` tuples. To retrieve data in the same
        order as it should be presented, iterate over the :attr:`rows` and
        :attr:`cols` attributes.
        """
        if self._data is None:
            self._data = {
                (row, col): self._func(row, col)
                for row in self.rows
                for col in self.cols
                }
            self._func = None
        return self._data

    @property
    def data_iter(self):
        """
        Returns an iterator over :attr:`data` where each key, value combination
        is returned as a tuple of (row_key, col_key, value), and each row and
        column key is returned as a tuple, regardless of the number of row and
        column dimensions. Furthermore, items are returned in declared row then
        column order. This property is intended to make renderers simpler.
        """
        for row_tuple, row_key in zip(self.rows_iter, self.rows):
            for col_tuple, col_key in zip(self.cols_iter, self.cols):
                yield row_tuple, col_tuple, self.data[(row_key, col_key)]

    @property
    def as_ndarray(self):
        """
        Returns the content of the table as a `numpy`_ :class:`ndarray` with
        the shape ``(rows, cols)``. Rows and columns will be in the order given
        by the :attr:`rows` and :attr:`cols` attributes. Please note that row
        and column keys are *not* included in the resulting array (as ndarrays
        purposely do not support heterogeneous data types).

        .. warning::

            Accessing this property will implicitly import the numpy module.
            This is not done during module import to avoid creating an
            explicit dependency on numpy.

        .. _numpy: http://www.numpy.org/
        """
        import numpy as np
        return np.asarray(
            [
                [self.data[(row, col)] for col in self.cols]
                for row in self.rows
                ], dtype=np.float)

    @property
    def as_dataframe(self):
        """
        Returns the content of the table as a `pandas`_ :class:`DataFrame`. The
        :attr:`rows` and :attr:`cols` attributes will be included as the index
        and columns of the resulting DataFrame.

        .. warning::

            Accessing this property will implicitly import the pandas module.
            This is not done during module import to avoid creating an
            explicit dependency on pandas.

        .. _pandas: http://pandas.pydata.org/
        """
        import pandas as pd
        if self.row_dims == 1:
            rows_index = pd.Index(self.rows, name=self.rows_title)
        else:
            rows_index = pd.MultiIndex.from_tuples(self.rows, names=self.rows_title)
        if self.col_dims == 1:
            cols_index = pd.Index(self.cols, name=self.cols_title)
        else:
            cols_index = pd.MultiIndex.from_tuples(self.cols, names=self.cols_title)
        return pd.DataFrame(self.as_ndarray, index=rows_index, columns=cols_index)

    def __repr__(self):
        return '<Table name="%s">' % self.name

    def __str__(self):
        # Calculate maximum columns lengths. Firstly, row headers
        if self.row_dims > 1:
            max_col_lens = [
                max(len(str(row[i]).strip()) for row in self.rows)
                for i in range(self.row_dims)
                ]
        else:
            max_col_lens = [
                max(len(str(row).strip()) for row in self.rows)
                ]
        # then column headers and data...
        if self.col_dims > 1:
            max_col_lens.extend([
                max(
                    max(len(str(self.data[(row, col)]).strip()) for row in self.rows), # lengths of all col values
                    *(len(str(header).strip()) for header in col) # length of col header(s)
                    )
                for col in self.cols
                ])
        else:
            max_col_lens.extend([
                max(
                    max(len(str(self.data[(row, col)]).strip()) for row in self.rows),
                    len(str(col).strip()),
                    )
                for col in self.cols
                ])

        result = ''
        # Print the column headers
        if self.col_dims > 1:
            for i in range(self.col_dims):
                result += ' | '.join(
                    '%*s' % (max_col_len, str(col[i]))
                    for max_col_len, col in zip(
                        max_col_lens, (('',) * self.col_dims,) * self.row_dims + self.cols)
                    )
                result += '\n'
        else:
            result += ' | '.join(
                '%*s' % (max_col_len, str(col))
                for max_col_len, col in zip(
                    max_col_lens, ('',) * self.row_dims + self.cols)
                )
            result += '\n'

        # Print a separator row
        result += '-+-'.join(
            '-' * max_col_len
            for max_col_len in max_col_lens
            )
        result += '\n'

        # Print the data rows
        for row in self.rows:
            if self.row_dims > 1:
                result += ' | '.join(
                    '%*s' % (max_col_len, value)
                    for max_col_len, value in zip(
                        max_col_lens,
                        [
                            str(head).strip()
                            for head in row
                            ] +
                        [
                            str(self.data[(row, col)]).strip()
                            for col in self.cols
                            ]
                        )
                    )
            else:
                result += ' | '.join(
                    '%*s' % (max_col_len, value)
                    for max_col_len, value in zip(
                        max_col_lens,
                        [str(row).strip()] +
                        [
                            str(self.data[(row, col)]).strip()
                            for col in self.cols
                            ]
                        )
                    )
            result += '\n'
        if sys.version_info.major < 3:
            return result.encode('utf-8')
        return result

