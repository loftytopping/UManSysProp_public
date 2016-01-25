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
An online framework for calculating the properties of individual organic
molecules and ensemble mixtures
"""

from __future__ import (
    unicode_literals,
    absolute_import,
    print_function,
    division,
    )
str = type('')

import pkgutil
import json
from textwrap import dedent

from flask import (
    Flask,
    request,
    url_for,
    render_template,
    make_response,
    send_file,
    jsonify,
    abort,
    )
import docutils.core

from . import tools
from . import renderers
from . import forms

app = Flask(__name__)
# maximum file upload is 1Mb
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024
# equalto was only added in Jinja 2.8 ?!
app.jinja_env.tests.setdefault('equalto', lambda value, other: value == other)

tools = {
    modname.split('.')[-1]: finder.find_module(modname).load_module(modname)
    for (finder, modname, ispkg) in pkgutil.iter_modules(
        tools.__path__, prefix=tools.__name__ + '.')
    if not ispkg and modname != 'template'
    }


@app.route('/')
@app.route('/<name>')
def template(name='index'):
    return render_template(
        '%s.html' % name,
        tools=tools,
        )


@app.route('/api')
def api():
    mimetype = request.accept_mimetypes.best_match([
        'text/html',
        'application/json',
        ])
    if mimetype == 'text/html':
        return template('api')
    elif mimetype == 'application/json':
        response = jsonify(**{
            mod_name: {
                'url': url_for('call', name=mod_name),
                'title': (mod.__doc__ or '').strip(),
                'doc': dedent(mod.handler.__doc__ or ''),
                'params': [
                    field.name for field in mod.HandlerForm()
                    if field.name not in ('csrf_token', 'output_format')
                    ],
                }
            for mod_name, mod in tools.items()
            })
        # Simple CORS setup
        response.headers['Access-Control-Allow-Origin'] = '*'
        return response
    else:
        abort(406)


@app.route('/api/<name>', methods=['GET'])
def api_docs(name):

    def render_docs(docstring):
        if not isinstance(docstring, str):
            docstring = docstring.decode('utf-8')
        docstring = dedent(docstring)
        result = docutils.core.publish_parts(
                docstring, writer_name='html', settings_overrides={
                    'input_encoding': 'unicode',
                    'output_encoding': 'unicode',
                    })
        return result['fragment']

    try:
        return render_template(
            'api_docs.html',
            title=name,
            name=name,
            tool=tools[name],
            render_docs=render_docs,
            )
    except KeyError:
        abort(404)


@app.route('/api/<name>', methods=['POST'])
def call(name):
    # Ensure CORS is on for all responses, including errors
    headers = {'Access-Control-Allow-Origin': '*'}
    # Fail if the RPC call has more than a meg of data
    if request.content_length > 1048576:
        result = jsonify(exc_type='ValueError', exc_value='Request too large')
        status = 413
    else:
        try:
            mod = tools[name]
        except KeyError:
            result = jsonify(exc_type='NameError', exc_value='Unknown method')
            status = 404
        else:
            try:
                args = json.loads(request.get_data(cache=False, as_text=True))
                args = forms.convert_args(mod.HandlerForm(formdata=None), args)
            except ValueError as e:
                result = jsonify(exc_type='ValueError', exc_value='Badly formed parameters: %s' % str(e))
                status = 400
            except KeyError as e:
                result = jsonify(exc_type='KeyError', exc_value='Missing parameter: %s' % str(e))
                status = 400
            else:
                try:
                    result = mod.handler(**args)
                except (ValueError, KeyError) as e:
                    result = jsonify(exc_type=e.__class__.__name__, exc_value=str(e))
                    status = 400
                else:
                    headers, result = renderers.render('application/json', result)
                    status = 200
    response = make_response(result)
    response.mimetype = 'application/json'
    response.headers.extend(headers)
    return response, status


@app.route('/tool/<name>', methods=['GET', 'POST'])
def tool(name):
    # Present the tool's input form, or execute the tool's handler callable
    # based on whether the HTTP request is a GET or a POST
    try:
        mod = tools[name]
    except KeyError:
        abort(404)
    form = mod.HandlerForm(request.form)
    if form.validate_on_submit():
        args = form.data
        mimetype = args.pop('output_format')
        result = mod.handler(**args)
        headers, result = renderers.render(mimetype, result)
        # If we're generating HTML, wrap the result in a template
        if mimetype == 'text/html':
            result = render_template(
                'result.html',
                title=mod.__doc__,
                result=result)
        response = make_response(result)
        response.mimetype = mimetype
        response.headers.extend(headers)
        return response
    return render_template(
        '%s.html' % name,
        title=mod.__doc__,
        form=form,
        )


def main():
    app.secret_key = 'testing'
    app.run(
        host='0.0.0.0',
        debug=True
        )
