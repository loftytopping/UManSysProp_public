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
    range = xrange
except NameError:
    pass
try:
    import copy_reg as copyreg
except ImportError:
    import copyreg


import math
import decimal
import itertools
from collections import namedtuple

import pybel
from flask import request
from flask_wtf import Form as DeclarativeForm
from wtforms.fields import (
    Field,
    BooleanField,
    FloatField as _FloatField,
    RadioField,
    SelectField,
    SelectMultipleField,
    StringField,
    PasswordField,
    TextAreaField,
    HiddenField,
    FileField,
    FormField,
    SubmitField,
    FieldList,
    )
from wtforms.fields.html5 import (
    SearchField,
    TelField,
    URLField,
    EmailField,
    DateField,
    DateTimeField,
    IntegerField,
    DecimalField,
    )
from wtforms.validators import (
    ValidationError,
    Optional,
    DataRequired,
    Email,
    EqualTo,
    IPAddress,
    MacAddress,
    Length,
    InputRequired,
    NumberRange,
    Regexp,
    URL,
    UUID,
    AnyOf,
    NoneOf,
    )
from wtforms.widgets import TextInput
from wtforms.widgets.html5 import NumberInput

from .html import html, literal, content, tag
from . import data
from . import groups
from . import renderers


class Form(DeclarativeForm):
    output_format = SelectField(
        'Output format', choices=renderers.registered(), default='text/html')


class SubForm(DeclarativeForm):
    def __init__(self, csrf_enabled=False, *args, **kwargs):
        super(SubForm, self).__init__(*args, csrf_enabled=False, **kwargs)


def unpickle_molecule(s):
    return pybel.readstring(b'smi', s)
def pickle_molecule(m):
    return unpickle_molecule, (str(m).strip().encode('ascii'),)
copyreg.pickle(pybel.Molecule, pickle_molecule)


def smiles(s):
    """
    Converts *s* into an OpenBabel :class:`Molecule` object. Raises
    :exc:`ValueError` if *s* is not a valid SMILES string.
    """
    if isinstance(s, str):
        s = s.encode('ascii')
    try:
        return pybel.readstring(b'smi', s)
    except IOError:
        raise ValueError('"%s" is not a valid SMILES string' % s)


def frange(start, stop=None, step=1.0):
    """
    Floating point variant of :func:`range`. Note that this variant has several
    inefficiencies compared to the built-in range, notably that reversal of
    the resulting generator relies enumeration of the generator.
    """
    if stop is None:
        stop, start = start, 0.0
    count = int(math.ceil((stop - start) / step))
    return (start + n * step for n in range(count))


class ZeroIonCharge(object):
    def __init__(self, message=None):
        if message is None:
            message = 'Non-zero charge balance'
        self.message = message

    def __call__(self, form, field):
        if groups.aggregate_matches(
            groups.all_matches(data.AIOMFAC_ION_SMARTS, field.data),
            data.AIOMFAC_ION_CHARGE
            ) != 0.0:
            raise ValidationError(self.message)


class FloatField(_FloatField):
    widget = NumberInput(step='any')


class SMILESField(Field):
    """
    Represents a text input which accepts a SMILES strings representing
    a chemical compound. The field's data is returned as an OpenBabel molecule
    object.

    :param compounds:
        If provided, a sequence of ``(value, label)`` tuples which can be
        selected by drop-down from the text field. Defaults to an empty
        sequence.
    """

    widget = TextInput()
    compounds = ()

    def __init__(self, label=None, validators=None, compounds=None, **kwargs):
        super(SMILESField, self).__init__(label, validators, **kwargs)
        if compounds is not None:
            self.compounds = compounds

    def __call__(self, **kwargs):
        if self.compounds:
            kwargs['list'] = '%s-list' % self.id
        return super(SMILESField, self).__call__(**kwargs)

    @property
    def scripts(self):
        return tag.datalist(
            (
                tag.option(value, label=label)
                for (value, label) in self.compounds
                ),
            id='%s-list' % self.id
            )

    def _value(self):
        if self.data:
            return self.data.write(b'smi').decode('ascii')
        else:
            return u''

    def process_formdata(self, valuelist):
        if valuelist:
            self.data = smiles(valuelist[0])
        else:
            self.data = None


class SMILESListField(FormField):
    """
    Represents a compound input which defines a list of SMILES strings
    representing chemical compounds, or accepts a file upload containing one
    SMILES string per line. The field's data is returned as a sequence of
    OpenBabel Molecule objects.

    Additional keyword arguments introduced by this class are:

    :param entry_label:
        Provides the label appearing above the SMILES text entry field

    :param upload_label:
        Provides the label appearing beside the file upload field

    :param compounds:
        If provided, a sequence of ``(value, label)`` tuples which can be
        selected by drop-down from the text-entry field. Defaults to an empty
        sequence.
    """

    def __init__(self, label=None, validators=None, entry_label=None,
            upload_label=None, compounds=None, **kwargs):

        if entry_label is None:
            entry_label = 'SMILES'
        if upload_label is None:
            upload_label = 'Upload SMILES'

        # The ZeroIonCharge validator applies to the whole list
        validators = validators or []
        self.validators = [
                v for v in validators
                if isinstance(v, (Length, ZeroIonCharge))]
        validators = [
                v for v in validators
                if not isinstance(v, (Length, ZeroIonCharge))]

        class ListForm(SubForm):
            entry = FieldList(
                SMILESField(entry_label, validators=validators, compounds=compounds))
            upload = FileField(upload_label)

        super(SMILESListField, self).__init__(ListForm, label, **kwargs)

    def validate(self, form, extra_validators=tuple()):
        self.form.validate()
        chain = itertools.chain(self.validators, extra_validators)
        self._run_validation_chain(form, chain)
        return len(self.errors) == 0

    @property
    def data(self):
        if self.form.upload.name in request.files:
            try:
                return [
                    smiles(line)
                    for i, _line in enumerate(
                        request.files[self.form.upload.name].read().splitlines(),
                        start=1
                        )
                    for line in (_line.strip(),)
                    if line and not line.startswith('#')
                    ]
            except ValueError as e:
                raise ValueError('%s on line %d' % (str(e), i))
        else:
            return self.form.entry.data

    def __call__(self, **kwargs):
        if not len(self.form.entry):
            self.form.entry.append_entry()
        # XXX Layout specific to UManSysProp
        return tag.div(
            tag.div(
                tag.div(
                    tag.label(
                        tag.input(
                            type='checkbox',
                            value='file'
                            ),
                        ' upload file',
                        ),
                    class_='medium-10 small-9 columns'
                    ),
                tag.div(
                    tag.a('Add', class_='button radius tiny right', data_toggle='fieldset-add-row'),
                    class_='medium-2 small-3 columns clearfix'
                    ),
                class_='row'
                ),
            tag.div(
                tag.div(
                    self.form.upload,
                    class_='small-12 columns'
                    ),
                class_='row',
                data_toggle='fieldset-upload'
                ),
            (tag.div(
                tag.div(
                    field,
                    class_='medium-10 small-9 columns'
                    ),
                tag.div(
                    tag.a('Remove', class_='button radius tiny right', data_toggle='fieldset-remove-row'),
                    class_='medium-2 small-3 columns clearfix'
                    ),
                class_='row',
                data_toggle='fieldset-entry'
                ) for field in self.form.entry),
            id=self.id,
            data_toggle='fieldset',
            data_freeid=len(self.form.entry)
            )

    @property
    def scripts(self):
        template = """\
$('div#%(id)s').each(function() {
    var $field = $(this);
    var $check = $field.find(':checkbox');
    var $add = $field.find('a[data-toggle=fieldset-add-row]');
    var $remove = $field.find('a[data-toggle=fieldset-remove-row]');
    var $upload = $field.find('div[data-toggle=fieldset-upload]');
    var freeid = parseInt($field.data('freeid'));

    $add.click(function() {
        // Find the last row and clone it
        var $oldrow = $field.find('div[data-toggle=fieldset-entry]:last');
        var $row = $oldrow.clone(true);
        // Re-write the ids of the input in the row
        $row.find(':input').each(function() {
            var newid = $(this).attr('id').replace(
                /%(id)s-entry-(\d{1,4})/,
                '%(id)s-entry-' + freeid);
            $(this)
                .attr('name', newid)
                .attr('id', newid)
                .val('')
                .removeAttr('checked');
        });
        $oldrow.after($row);
        freeid++;
    });

    $remove.click(function() {
        if ($field.find('div[data-toggle=fieldset-entry]').length > 1) {
            var thisRow = $(this).closest('div[data-toggle=fieldset-entry]');
            thisRow.remove();
        }
    });

    $check.change(function() {
        // Refresh the entry matches
        var $entry = $field.find('div[data-toggle=fieldset-entry]').add($add);

        var $show = this.checked ? $upload : $entry;
        var $hide = this.checked ? $entry : $upload;
        $hide.fadeOut('fast', function() {
            $hide
                .find(':input')
                    .prop('disabled', true);
            $show
                .find(':input')
                    .prop('disabled', false)
                .end()
                    .fadeIn('fast');
        });
    });

    if ($field.find(':file').val()) {
        $check.prop('checked', true);
        $field.find('div[data-toggle=fieldset-entry]').add($add)
            .hide().find(':input').prop('disabled', true);
    }
    else {
        $check.prop('checked', false);
        $upload
            .hide().find(':input').prop('disabled', true);
    }
});
"""
        return literal('\n'.join(
            [tag.script(literal(template % {'id': self.id}))] +
            [field.scripts for field in self.form.entry]
            ))


class SMILESDictField(FormField):
    """
    Represents a compound input which defines a mapping of SMILES strings to
    floating point values, or accepts a file upload containing one SMILES
    string and one floating point value separated by whitespace, per line. The
    field's data is returned as a mapping of OpenBabel Molecule objects to
    float values.

    Additional keyword arguments introduced by this class are:

    :param entry_label:
        Provides the label appearing above the SMILES text entry field

    :param data_label:
        Provides the label appearing above the floating point entry field

    :param upload_label:
        Provides the label appearing beside the file upload field

    :param compounds:
        If provided, a sequence of ``(value, label)`` tuples which can be
        selected by drop-down from the text-entry field. Defaults to an empty
        sequence.
    """

    def __init__(self, label=None, validators=None, entry_label=None,
            data_label=None, upload_label=None, compounds=None, **kwargs):

        if entry_label is None:
            entry_label = 'SMILES'
        if data_label is None:
            data_label = 'Concentration'
        if upload_label is None:
            upload_label = 'Upload SMILES'

        # The Length and ZeroIonCharge validators applies to the whole dict
        validators = validators or []
        self.validators = [
                v for v in validators
                if isinstance(v, (Length, ZeroIonCharge))]
        validators = [
                v for v in validators
                if not isinstance(v, (Length, ZeroIonCharge))]

        # If a NumberRange validator has been specified, extract its minimums
        # and maximums
        self._min = None
        self._max = None
        for v in validators:
            if isinstance(v, NumberRange):
                self._min = v.min
                self._max = v.max
                break

        class MapForm(SubForm):
            smiles = SMILESField(entry_label,
                    validators=[
                        v for v in validators
                        if not isinstance(v, NumberRange)],
                    compounds=compounds)
            data = FloatField(data_label,
                    validators=[
                        v for v in validators
                        if isinstance(v, NumberRange)])

        class ListForm(SubForm):
            entry = FieldList(FormField(MapForm))
            upload = FileField(upload_label)

        super(SMILESDictField, self).__init__(ListForm, label, **kwargs)

    def validate(self, form, extra_validators=tuple()):
        self.form.validate()
        chain = itertools.chain(self.validators, extra_validators)
        self._run_validation_chain(form, chain)
        return len(self.errors) == 0

    @property
    def data(self):
        if self.form.upload.name in request.files:
            result = {}
            try:
                # XXX Check each associated value is >=0
                for i, line in enumerate(
                        request.files[self.form.upload.name].read().splitlines(),
                        start=1):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        key, value = line.split(None, 1)
                        result[smiles(key)] = float(value)
                return result
            except ValueError as e:
                e.args += ('on line %d' % i)
                raise
        else:
            return {
                e.smiles.data: e.data.data
                for e in self.form.entry
                }

    def __call__(self, **kwargs):
        if not len(self.form.entry):
            self.form.entry.append_entry()
        # XXX Layout specific to UManSysProp
        return tag.div(
            tag.div(
                tag.div(
                    tag.label(
                        tag.input(
                            type='checkbox',
                            value='file'
                            ),
                        ' upload file',
                        ),
                    class_='small-12 columns'
                    ),
                class_='row'
                ),
            tag.div(
                tag.div(
                    self.form.upload,
                    class_='small-12 columns'
                    ),
                class_='row',
                data_toggle='fieldset-upload'
                ),
            tag.div(
                tag.div(
                    self.form.entry[0].smiles.label(class_='inline'),
                    class_='small-6 columns'
                    ),
                tag.div(
                    self.form.entry[0].data.label(class_='inline', min=self._min, max=self._max),
                    class_='medium-4 small-3 columns'
                    ),
                tag.div(
                    tag.a('Add', class_='button radius tiny right', data_toggle='fieldset-add-row'),
                    class_='medium-2 small-3 columns clearfix'
                    ),
                class_='row',
                data_toggle='fieldset-entry'
                ),
            (tag.div(
                tag.div(
                    entry.smiles,
                    class_='small-6 columns'
                    ),
                tag.div(
                    entry.data,
                    class_='medium-4 small-3 columns'
                    ),
                tag.div(
                    tag.a('Remove', class_='button radius tiny right', data_toggle='fieldset-remove-row'),
                    class_='medium-2 small-3 columns clearfix'
                    ),
                class_='row',
                data_toggle='fieldset-entry'
                ) for entry in self.form.entry),
            id=self.id,
            data_toggle='fieldset',
            data_freeid=len(self.form.entry)
            )

    @property
    def scripts(self):
        template = """\
$('div#%(id)s').each(function() {
    var $field = $(this);
    var $check = $field.find(':checkbox');
    var $add = $field.find('a[data-toggle=fieldset-add-row]');
    var $remove = $field.find('a[data-toggle=fieldset-remove-row]');
    var $upload = $field.find('div[data-toggle=fieldset-upload]');
    var freeid = parseInt($field.data('freeid'));

    $add.click(function() {
        // Find the last row and clone it
        var $oldrow = $field.find('div[data-toggle=fieldset-entry]:last');
        var $row = $oldrow.clone(true);
        // Re-write the ids of the input in the row
        $row.find(':input').each(function() {
            var newid = $(this).attr('id').replace(
                /%(id)s-entry-(\d{1,4})/,
                '%(id)s-entry-' + freeid);
            $(this)
                .attr('name', newid)
                .attr('id', newid)
                .val('')
                .removeAttr('checked');
        });
        $oldrow.after($row);
        freeid++;
    });

    $remove.click(function() {
        if ($field.find('div[data-toggle=fieldset-entry]').length > 2) {
            var thisRow = $(this).closest('div[data-toggle=fieldset-entry]');
            thisRow.remove();
        }
    });

    $check.change(function() {
        // Refresh the entry matches
        var $entry = $field.find('div[data-toggle=fieldset-entry]');

        var $show = this.checked ? $upload : $entry;
        var $hide = this.checked ? $entry : $upload;
        $hide.fadeOut('fast', function() {
            $hide
                .find(':input')
                    .prop('disabled', true);
            $show
                .find(':input')
                    .prop('disabled', false)
                .end()
                    .fadeIn('fast');
        });
    });

    if ($field.find(':file').val()) {
        $check.prop('checked', true);
        $field.find('div[data-toggle=fieldset-entry]')
            .hide().find(':input').prop('disabled', true);
    }
    else {
        $check.prop('checked', false);
        $upload
            .hide().find(':input').prop('disabled', true);
    }
});
"""
        return literal('\n'.join(
            [tag.script(literal(template % {'id': self.id}))] +
            [entry.smiles.scripts for entry in self.form.entry]
            ))


class FloatRangeField(FormField):
    """
    Represents a complex input which defines either a single floating point
    value, or a range of floating point values evenly spaced between two
    inclusive end-points. In either case, the field returns a sequence of
    floating point values as its data.
    """

    def __init__(self, label=None, validators=None, **kwargs):

        # The Length validator applies to the whole list. If one is present
        # we extract its minimums and maximums to use as the min and max
        # attributes of the count field
        validators = validators or []
        self.validators = [
                v for v in validators
                if isinstance(v, Length)]
        validators = [
                v for v in validators
                if not isinstance(v, Length)]
        self._min_len = None
        self._max_len = None
        if self.validators:
            self._min_len = self.validators[0].min
            self._max_len = self.validators[0].max

        # If a NumberRange validator has been specified, extract its minimums
        # and maximums
        self._min = None
        self._max = None
        for v in validators:
            if isinstance(v, NumberRange):
                self._min = v.min
                self._max = v.max
                break

        class RangeForm(SubForm):
            count = IntegerField(
                '', default=1,
                validators=[NumberRange(min=self._min_len, max=self._max_len)]
                )
            start = FloatField(
                'values from', default=kwargs.get('default'),
                validators=validators)
            stop = FloatField(
                'to', default=kwargs.get('default'),
                validators=validators)

            def validate(self):
                if not super(RangeForm, self).validate():
                    return False
                if self.stop.data is not None and self.start.data > self.stop.data:
                    self.start.errors.append(
                        'Starting value must be less than ending value')
                    return False
                return True

        super(FloatRangeField, self).__init__(RangeForm, label, validators=None, **kwargs)

    def validate(self, form, extra_validators=tuple()):
        self.form.validate()
        chain = itertools.chain(self.validators, extra_validators)
        self._run_validation_chain(form, chain)
        return len(self.errors) == 0

    def __call__(self, **kwargs):
        # XXX Eliminate the range checkbox if self._min_len > 1
        return tag.div(
            tag.div(
                tag.div(
                    tag.label(
                        tag.input(
                            type='checkbox',
                            value='range'
                            ),
                        ' range'
                        ),
                    class_='small-12 columns'
                    ),
                class_='row'
                ),
            tag.div(
                tag.div(
                    self.form.start(id=self.form.start.id + '-single', min=self._min, max=self._max),
                    self.form.stop(id=self.form.stop.id + '-single', type='hidden'),
                    self.form.count(id=self.form.count.id + '-single', type='hidden', value=1),
                    class_='small-12 columns'
                    ),
                class_='row',
                data_toggle='single'
                ),
            tag.div(
                tag.div(
                    self.form.count(min=self._min_len, max=self._max_len),
                    class_='medium-3 columns',
                    data_toggle='count'
                    ),
                tag.div(
                    self.form.start.label(class_='inline'),
                    class_='medium-2 columns medium-text-center'
                    ),
                tag.div(
                    self.form.start(min=self._min, max=self._max),
                    class_='medium-3 columns'
                    ),
                tag.div(
                    self.form.stop.label(class_='inline'),
                    class_='medium-1 columns medium-text-center'
                    ),
                tag.div(
                    self.form.stop(min=self._min, max=self._max),
                    class_='medium-3 columns'
                    ),
                class_='row',
                data_toggle='multi'
                ),
            id=self.id
            )

    @property
    def scripts(self):
        template = """\
$('div#%(id)s').each(function() {
    var $field = $(this);
    var $check = $field.find(':checkbox');
    var $single = $field.find('div[data-toggle=single]');
    var $multi = $field.find('div[data-toggle=multi]');

    $single.find('#%(id)s-start-single').change(function() {
        $single.find('#%(id)s-stop-single').val($(this).val());
    });

    $check.change(function() {
        $show = this.checked ? $multi : $single;
        $hide = this.checked ? $single : $multi;
        $hide.fadeOut('fast', function() {
            $hide
                .find(':input')
                    .prop('disabled', true);
            $show
                .find(':input')
                    .prop('disabled', false)
                .end()
                    .fadeIn('fast');
        });
    });

    if (parseInt($multi.find('#%(id)s-count').val()) > 1) {
        $check.prop('checked', true);
        $single.hide().find(':input').prop('disabled', true);
        $multi.find(':input').prop('disabled', false);
    }
    else {
        $check.prop('checked', false);
        $multi.hide().find(':input').prop('disabled', true);
        $single.find(':input').prop('disabled', false);
    }

});
"""
        return tag.script(literal(template % {'id': self.id}))

    @property
    def data(self):
        start = self.form.start.data
        stop = self.form.stop.data
        count = self.form.count.data
        if count == 1:
            return [start]
        else:
            step = (stop - start) / (count - 1)
            return list(frange(start, stop, step)) + [stop]


class CoreAbundance(namedtuple('CoreAbundance', (
    'amount',
    'weight',
    'dissociation',
    ))):

    @property
    def concentration(self):
        return (self.amount / self.weight) * self.dissociation


class CoreAbundanceField(FormField):
    """
    Represents a complex input which defines an inert core abundance and its
    molecular weight. If soluble is True (which it is by default), the extra
    dissociation factor field is included in the interface. If it is False,
    the field is excluded from the interface, and a fixed value of 1.0 is
    returned in the namedtuple provided by the data property.
    """

    def __init__(self, label=None, validators=None, soluble=True, **kwargs):
        validators = validators or []
        self.soluble = soluble

        class AbundanceForm(SubForm):
            amount = FloatField(
                'Amount (µg/m³)', default=0.0,
                validators=[NumberRange(min=0.0)] + validators
                )
            weight = FloatField(
                'Weight (g/mol)', default=320.0,
                validators=[NumberRange(min=0.0)] + validators
                )
            dissociation = FloatField(
                'Dissociation factor (unitless)', default=1.0,
                validators=[NumberRange(min=0.0)] + validators
                )

        if not self.soluble:
            del AbundanceForm.dissociation
        super(CoreAbundanceField, self).__init__(
            AbundanceForm, label, validators=None, **kwargs)

    def __call__(self, **kwargs):
        return tag.div(
            tag.div(
                tag.div(
                    self.form.amount.label(class_='inline'),
                    class_='medium-2 columns',
                    ),
                tag.div(
                    self.form.amount(min=0.0),
                    class_='medium-2 columns',
                    ),
                tag.div(
                    self.form.weight.label(class_='inline'),
                    class_='medium-2 columns',
                    ),
                tag.div(
                    self.form.weight(min=0.0),
                    class_='medium-2 columns',
                    ),
                tag.div(
                    self.form.dissociation.label(class_='inline'),
                    class_='medium-2 columns',
                    ) if self.soluble else tag.div(class_='medium-4 columns'),
                tag.div(
                    self.form.dissociation(min=0.0),
                    class_='medium-2 columns',
                    ) if self.soluble else '',
                class_='row'
                ),
            id=self.id
            )

    @property
    def data(self):
        return CoreAbundance(
            amount=self.form.amount.data,
            weight=self.form.weight.data,
            dissociation=self.form.dissociation.data if self.soluble else 1.0,
            )

class Sizedependance(namedtuple('Sizedependance', (
    'diameter',
    'surface_tension',
    ))):

    @property
    def size(self):
        return self.diameter


class SizedependanceField(FormField):
    """
    Represents a complex input which defines a dry size and surface tension.
   
    """

    def __init__(self, label=None, validators=None, **kwargs):
        validators = validators or []

        class SizeForm(SubForm):
            diameter = FloatField(
                'Dry diameter (nm)', default=100.0,
                validators=[NumberRange(min=10.0)] + validators
                )
            surface_tension = FloatField(
                'Surface tension (mN/m)', default=72.0,
                validators=[NumberRange(min=10.0,max=140.0)] + validators
                )

        super(SizedependanceField, self).__init__(
            SizeForm, label, validators=None, **kwargs)

    def __call__(self, **kwargs):
        return tag.div(
            tag.div(
                tag.div(
                    self.form.diameter.label(class_='inline'),
                    class_='medium-2 columns',
                    ),
                tag.div(
                    self.form.diameter(min=10.0),
                    class_='medium-2 columns',
                    ),
                tag.div(
                    self.form.surface_tension.label(class_='inline'),
                    class_='medium-2 columns',
                    ),
                tag.div(
                    self.form.surface_tension(min=20.0,max=140.0),
                    class_='medium-2 columns',
                    ),
              class_='row'
                ),
            id=self.id
            )

    @property
    def data(self):
        return Sizedependance(
            diameter=self.form.diameter.data,
            surface_tension=self.form.weight.data,
            )


def convert_args(form, args):
    """
    Given a *form* and a dictionary of *args* which has been decoded from JSON,
    returns *args* with the type of each value converted for the corresponding
    field.
    """
    conversion = {
        IntegerField:       int,
        FloatField:         float,
        DecimalField:       decimal.Decimal,
        BooleanField:       bool,
        SMILESField:        smiles,
        SMILESListField:    lambda l: [smiles(s) for s in l],
        SMILESDictField:    lambda l: {smiles(s): float(v) for (s, v) in l.items()},
        FloatRangeField:    lambda l: [float(f) for f in l],
        CoreAbundanceField: lambda l: CoreAbundance(*l),
        }
    return {
        field.name: conversion.get(field.__class__, lambda x: x)(args[field.name])
        for field in form
        if field.name not in ('csrf_token', 'output_format')
        }

