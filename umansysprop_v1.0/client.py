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
=============================
``umansysprop.client`` Module
=============================

This module contains the client library for interacting with the UManSysProp
server. Only one user-accessible class is defined in the module:

UManSysProp
===========

.. autoclass:: UManSysProp

"""

from __future__ import (
    unicode_literals,
    absolute_import,
    print_function,
    division,
    )
str = type('')


import types
import json
try:
    from urllib.parse import urljoin
except ImportError:
    from urlparse import urljoin

import requests

from . import results


class UManSysProp(object):
    """
    Provides a simple Python interface to the methods provided via the `JSON
    API`_ on the `UManSysProp`_ website. Constructing an instance of this class
    will cause the new instance to query the server for all available methods.
    Each method will be exposed as a method of the instance, with docstrings
    obtained from the server. For example::

        >>> import umansysprop.client
        >>> client = umansysprop.client.UManSysProp()
        >>> client.vapour_pressure
        <bound method ?.vapour_pressure of <umansysprop.client.UManSysProp object at 0x7f80fbae7b50>>
        >>> client.vapour_pressure.__doc__
        u"\\nCalculates vapour pressures for all specified *compounds* (given..."

    The class can be constructed with an alternative *base_url* if you wish to
    point it a different server. The *base_url* parameter defaults to the
    `UManSysProp`_ website.

    Methods can be called like a normal Python method, but will result in a
    request being sent to the web-server, processed, and the JSON-formatted
    results being re-constructed as a :class:`~umansysprop.results.Result` list
    on the client.  For example::

        >>> result = client.vapour_pressure([
        ... 'CCCC', 'C(CC(=O)O)C(=O)O', 'C(=O)(C(=O)O)O',
        ... 'CCCCC/C=C/C/C=C/CC/C=C/CCCC(=O)O'],
        ... [298.15, 299.15, 300.15, 310.15],
        ... 'nannoolal', 'nannoolal')
        >>> result
        [<Table name="pressures">]
        >>> result.pressures
        <Table name="pressures">
        >>> result.pressures.data
        {(298.15, u'C(=O)(C(=O)O)O'): -5.196360545314141,
         (298.15, u'C(CC(=O)O)C(=O)O'): -6.33293991047814,
         (298.15, u'CCCC'): 0.22091492301164387,
         (298.15, u'CCCCC/C=C/C/C=C/CC/C=C/CCCC(=O)O'): -9.660331395164858,
         (299.15, u'C(=O)(C(=O)O)O'): -5.151703772558254,
         (299.15, u'C(CC(=O)O)C(=O)O'): -6.281177618554678,
         (299.15, u'CCCC'): 0.2354793193478599,
         (299.15, u'CCCCC/C=C/C/C=C/CC/C=C/CCCC(=O)O'): -9.58901500825095,
         (300.15, u'C(=O)(C(=O)O)O'): -5.107428775107059,
         (300.15, u'C(CC(=O)O)C(=O)O'): -6.229864995169396,
         (300.15, u'CCCC'): 0.24993365754879804,
         (300.15, u'CCCCC/C=C/C/C=C/CC/C=C/CCCC(=O)O'): -9.518352766693454,
         (310.15, u'C(=O)(C(=O)O)O'): -4.684643528880829,
         (310.15, u'C(CC(=O)O)C(=O)O'): -5.74023509658808,
         (310.15, u'CCCC'): 0.38868830156274703,
         (310.15, u'CCCCC/C=C/C/C=C/CC/C=C/CCCC(=O)O'): -8.845815136269506}

    Please refer to the reference for :class:`~umansysprop.results.Result` and
    :class:`~umansysprop.results.Table` for more information on accessing the
    result data.

    .. _UManSysProp: http://umansysprop.seaes.manchester.ac.uk/
    .. _JSON API: http://umansysprop.seaes.manchester.ac.uk/api
    """

    _method_template = """\
def {name}(self, {params}):
    return self._json_rpc("{url}", {call})
"""

    def __init__(self, base_url='http://umansysprop.seaes.manchester.ac.uk/'):
        self._base_url = base_url
        response = requests.get(urljoin(self._base_url, 'api'), headers={
            'Accept': 'application/json'})
        response.raise_for_status()
        for name, props in response.json().items():
            # Construct a dynamic method for each function that the API
            # defines.  The method will take the parameters specified by the
            # API, and will have a doc-string also specified by the API
            method_definition = self._method_template.format(
                    name=name,
                    url=props['url'],
                    params=', '.join(props['params']),
                    call=', '.join(
                        '%s=%s' % (param, param) for param in props['params'])
                    )
            l = {}
            exec(method_definition, globals(), l)
            f = l[name]
            f.__doc__ = props['doc']
            setattr(self, name, types.MethodType(f, self))

    def _json_rpc(self, url, **params):
        response = requests.post(
                urljoin(self._base_url, url),
                data=json.dumps(params),
                headers={
                    'Accept': 'application/json',
                    'Content-Type': 'application/json',
                    })
        if 400 <= response.status_code < 500:
            # Some kind of client error; try and decode the body as JSON to
            # determine the details and raise a reasonable exception
            exc_type = response.json()['exc_type']
            exc_value = response.json()['exc_value']
            # Only permit a specific set of exceptions
            exc_class = {
                'ValueError':   ValueError,
                'NameError':    NameError,
                'KeyError':     KeyError,
                }[exc_type]
            if isinstance(exc_value, str):
                raise exc_class(exc_value)
            else:
                raise exc_class(*exc_value)
        elif response.status_code >= 500:
            raise RuntimeError('Server error: %s' % response.body)
        else:
            response.raise_for_status()
        return results.Result.from_json(response.json())

