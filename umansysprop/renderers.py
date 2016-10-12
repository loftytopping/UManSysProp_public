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

from __future__ import (
    unicode_literals,
    absolute_import,
    print_function,
    division,
    )
str = type('')
try:
    from itertools import izip as zip
except ImportError:
    pass


import sys
import io
import csv
import pickle
import tempfile
import textwrap

import xlsxwriter as xl
import pybel
from flask import json

from .zip import ZipFile, ZIP_DEFLATED
from .html import TagFactory


_RENDERERS = {}

def register(mimetype, label, headers=None):
    if headers is None:
        headers = {}
    def decorator(func):
        if mimetype in _RENDERERS:
            raise ValueError('A handler for MIME-type %s already exists' % mimetype)
        _RENDERERS[mimetype] = (label, headers, func)
        return func
    return decorator

def registered():
    return [
        (mimetype, label)
        for (mimetype, (label, headers, func)) in _RENDERERS.items()
        ]

def render(mimetype, obj, **kwargs):
    try:
        label, headers, func = _RENDERERS[mimetype]
    except KeyError:
        raise ValueError('Unknown MIME-type %s' % mimetype)
    else:
        return headers, func(obj, **kwargs)


def _format_key(value):
    # This rather hacky routine is here to deal with the crappy string
    # conversion from OpenBabel's Molecule class
    if isinstance(value, tuple):
        return tuple(_format_key(key) for key in value)
    elif isinstance(value, pybel.Molecule):
        return str(value).strip()
    else:
        return value


@register('application/json', 'JSON file', headers={
    # Simple CORS setup; see http://enable-cors.org/server.html for more
    # information
    'Access-Control-Allow-Origin': '*',
    })
def render_json(results, **kwargs):

    def render_table(table):
        return {
            'name': table.name,
            'title': table.title,
            'rows_title': table.rows_title,
            'cols_title': table.cols_title,
            'rows_unit': table.rows_unit,
            'cols_unit': table.cols_unit,
            'data': [
                {
                    'key': (_format_key(row_key), _format_key(col_key)),
                    'value': table.data[(row_key, col_key)],
                    }
                for row_key in table.rows
                for col_key in table.cols
                ]
            }

    return json.dumps([render_table(table) for table in results], **kwargs)


@register('application/octet-stream', 'Python pickle')
def render_pickle(results, **kwargs):
    for table in results:
        # Force evaluation of the data and wipe the function reference (it's
        # no longer needed and typically prevents serialization by being a
        # lambda of some sorts)
        table.data
    return pickle.dumps(results, protocol=0)


@register('application/xml', 'XML file', headers={
        'Content-Disposition': 'attachment; filename=umansysprop.xml',
        })
def render_xml(results, **kwargs):
    tag = TagFactory(xml=True)

    def render_table(table):
        return tag.table(
            tag.cols(
                tag.dim(title=title, unit=unit)
                for (title, unit) in zip(table.col_titles, table.col_units)
                ),
            tag.rows(
                tag.dim(title=title, unit=unit)
                for (title, unit) in zip(table.row_titles, table.row_units)
                ),
            tag.data(
                tag.datum(
                    tag.row(tag.dim(key=key) for key in row_key),
                    tag.col(tag.dim(key=key) for key in col_key),
                    value
                    )
                for row_key, col_key, value in table.data_iter
                ),
            title=table.title,
            name=table.name,
            )

    return tag.tables(render_table(table) for table in results)


@register('application/zip', 'Zipped CSV files', headers={
        'Content-Disposition': 'attachment; filename=umansysprop.zip',
        })
def render_csv(results, **kwargs):

    def render_table(table):
        # Deal with incompatibility between Py2 and Py3's CSV writer
        if sys.version_info.major == 3:
            stream = io.StringIO(newline='')
        else:
            stream = io.BytesIO()
        writer = csv.writer(stream)
        for dim in range(table.col_dims):
            writer.writerow(
                [''] * table.row_dims +
                [_format_key(col_key[dim]) for col_key in table.cols_iter]
                )
        for data_row, row_keys in zip(table.rows, table.rows_iter):
            writer.writerow(
                [_format_key(row_key) for row_key in row_keys] +
                [table.data[(data_row, data_col)] for data_col in table.cols]
                )
        stream.seek(0)
        return stream

    def render_readme():
        s = """\
This file details the CSV files contained within this archive, and the
structure of their rows and columns.
"""
        for table in results:
            s += '\n'
            s += '%s.csv:\n' % table.name
            s += ''.join(
                '  %s\n' % l for l in textwrap.wrap(table.title, width=70)
                )
            s += '  Rows:\n'
            for title, unit in zip(table.row_titles, table.row_units):
                if unit:
                    s += '    %s [%s]\n' % (title, unit if unit else 'unitless')
                else:
                    s += '    %s\n' % title
            s += '  Cols:\n'
            for title, unit in zip(table.col_titles, table.col_units):
                if unit:
                    s += '    %s [%s]\n' % (title, unit if unit else 'unitless')
                else:
                    s += '    %s\n' % title
        stream = io.BytesIO()
        stream.write(s.encode('utf-8'))
        stream.seek(0)
        return stream

    with io.BytesIO() as stream:
        with ZipFile(stream, 'w', compression=ZIP_DEFLATED) as archive:
            archive.comment = '\n'.join(table.title for table in results).encode('utf-8')
            for table in results:
                archive.write(render_table(table), '%s.csv' % table.name)
            archive.write(render_readme(), 'README.txt')
        return stream.getvalue()


@register('application/vnd.openxmlformats-officedocument.spreadsheetml.sheet', 'Excel file', headers={
        'Content-Disposition': 'attachment; filename=umansysprop.xlsx',
        })
def render_xlsx(results, **kwargs):
    stream = io.BytesIO()
    workbook = xl.Workbook(stream, {'in_memory': True})
    col_title_f = workbook.add_format({'bold': True, 'left': 1})
    row_title_f = workbook.add_format({'bold': True, 'bottom': 1})
    first_col_title_f = workbook.add_format({'bottom': 1, 'left': 1})
    first_data_f = workbook.add_format({'top': 1, 'left': 1})
    col_key_f = workbook.add_format({'top': 1})
    row_key_f = workbook.add_format({'left': 1})
    heading_f = workbook.add_format({
        'font_size': 14,
        'center_across': True,
        })

    def render_table(table):
        worksheet = workbook.add_worksheet(table.name[:31])
        # Write worksheet title
        worksheet.set_row(0, 24)
        worksheet.write(0, 0, table.title, heading_f)
        for i in range(table.row_dims + len(table.cols) - 1):
            worksheet.write_blank(0, i + 1, heading_f)
        # Write column and row titles
        column_titles = ' / '.join(
            '%s [%s]' % (title, unit) if unit else title
            for (title, unit) in zip(table.col_titles, table.col_units)
            )
        worksheet.merge_range(
            2, table.row_dims, 2, table.row_dims + len(table.cols) - 1,
            column_titles, col_title_f)
        for col_ix, (title, unit) in enumerate(zip(table.row_titles, table.row_units)):
            worksheet.write(
                table.col_dims + 2, col_ix, '%s [%s]' % (title, unit) if unit else title,
                row_title_f)
        # Write column keys
        for col_ix, col_key in enumerate(table.cols_iter, start=table.row_dims):
            span = table.col_spans[col_key]
            for col_dim in range(table.col_dims):
                if span[col_dim] > 1:
                    worksheet.merge_range(
                        col_dim + 3, col_ix,
                        col_dim + 3, col_ix + span[col_dim] - 1,
                        _format_key(col_key[col_dim]),
                        first_col_title_f if col_ix == table.row_dims else None)
                elif span[col_dim] == 1:
                    worksheet.write(
                        col_dim + 3, col_ix,
                        _format_key(col_key[col_dim]),
                        first_col_title_f if col_ix == table.row_dims else None)
        # Write row keys
        for row_ix, row_key in enumerate(table.rows_iter, start=table.col_dims + 3):
            span = table.row_spans[row_key]
            for row_dim in range(table.row_dims):
                if span[row_dim] > 1:
                    worksheet.merge_range(
                        row_ix, row_dim,
                        row_ix + span[row_dim] - 1, row_dim,
                        _format_key(row_key[row_dim]))
                elif span[row_dim] == 1:
                    worksheet.write(
                        row_ix, row_dim,
                        _format_key(row_key[row_dim]))
        # Write data
        for row_ix, row_key in enumerate(table.rows, start=table.col_dims + 3):
            for col_ix, col_key in enumerate(table.cols, start=table.row_dims):
                worksheet.write(
                    row_ix, col_ix, table.data[(row_key, col_key)],
                    first_data_f if (row_ix, col_ix) == (table.col_dims + 3, table.row_dims) else
                    col_key_f if row_ix == table.col_dims + 3 else
                    row_key_f if col_ix == table.row_dims else
                    None
                    )

    for table in results:
        render_table(table)
    workbook.close()
    return stream.getvalue()


@register('text/html', 'HTML (view in web browser)')
def render_html(results, **kwargs):
    tag = TagFactory(xml=False)

    def render_table(table):
        column_titles = ' / '.join(
            '%s [%s]' % (title, unit) if unit else title
            for (title, unit) in zip(table.col_titles, table.col_units)
            )
        return tag.table(
            tag.caption(table.title),
            tag.thead(
                tag.tr(
                    (tag.th('') for i in range(table.row_dims)),
                    tag.th(column_titles, colspan=len(table.cols)),
                    ),
                (
                    tag.tr(
                        (tag.th('') for i in range(table.row_dims))
                        if col_dim < table.col_dims - 1 else
                        (
                            tag.th('%s [%s]' % (row_title, row_unit) if row_unit else row_title)
                            for (row_title, row_unit) in zip(table.row_titles, table.row_units)
                            ),
                        (
                            tag.th(_format_key(key[col_dim]), colspan=span if span > 1 else None)
                            for key in table.cols_iter
                            for span in (table.col_spans[key][col_dim],)
                            if span > 0
                            )
                        )
                    for col_dim in range(table.col_dims)
                    )
                ),
            tag.tbody(
                tag.tr(
                    (
                        tag.th(_format_key(row_key[row_dim]), rowspan=span if span > 1 else None)
                        for row_dim in range(table.row_dims)
                        for span in (table.row_spans[row_key][row_dim],)
                        if span > 0
                        ),
                    (
                        tag.td(table.data[(data_row, data_col)])
                        for (data_col, col_key) in zip(table.cols, table.cols_iter)
                        )
                    )
                for (data_row, row_key) in zip(table.rows, table.rows_iter)
                ),
            id=table.name,
            )

    return tag.div(render_table(table) for table in results)

